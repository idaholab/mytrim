# MyTRIM Documentation

## Table of Contents
1. [Introduction](#introduction)
2. [Architecture and Design](#architecture-and-design)
3. [Using runmytrim](#using-runmytrim)
4. [Using runstopping](#using-runstopping)
5. [Customizing MyTRIM](#customizing-mytrim)
6. [Building and Installation](#building-and-installation)

---

## Introduction

MyTRIM is a three-dimensional binary collision Monte Carlo (BCA) library for simulating ion collision cascades in materials. It is based on the SRIM (Stopping and Range of Ions in Matter) methodology described in *"The Stopping and Range of Ions in Solids"* by J. F. Ziegler, J. P. Biersack, and U. Littmark.

Unlike the original TRIM/SRIM codes which were primarily 1D, MyTRIM is fully three-dimensional and supports arbitrary sample geometries including layered structures, wires, clusters, and time-dependent materials.

### Key Features
- Full 3D collision cascade simulation
- Multiple scattering potentials (Universal, Moliere, CKR)
- Extensible sample geometries
- Multi-threaded execution
- JSON-based input format
- Built-in stopping power calculations

---

## Architecture and Design

### Core Class Hierarchy

#### SimconfType (simconf.h:47)
The `SimconfType` class manages global simulation configuration and data tables.

**Key Responsibilities:**
- Loads ZBL coefficient tables from data files (SCOEF.95A/B, SNUC03.dat, ELNAME.dat)
- Manages random number generation using C++11 MT19937
- Tracks simulation statistics (vacancies created, energy losses)
- Handles length scale conversions

**Important Members:**
```cpp
Real ed;              // Displacement energy threshold
Real tmin;            // Minimum energy threshold
Real tau;             // Time constant
Real cw;              // Cascade width parameter
Real EelTotal;        // Total electronic energy loss
Real EnucTotal;       // Total nuclear energy loss
int vacancies_created; // Count of vacancies

Real drand();         // Random double [0,1)
unsigned int irand(); // Random integer
void seed(unsigned int seed);           // Set RNG seed
void setLengthScale(Real l);            // Set length units
```

**Length Scale:**
- `1` → Units are Angstroms (default)
- `10` → Units are Nanometers
- `10000` → Units are Micrometers

#### IonBase (ion.h:30)
Represents an ion or recoil atom in the simulation.

**Key Members:**
```cpp
int _Z;              // Atomic number
Real _m;             // Mass in amu
Real _E;             // Kinetic energy in eV
Point _dir;          // Normalized velocity vector
Point _pos;          // Position
unsigned int _seed;  // RNG seed for this ion
int _gen;            // Generation (0=PKA, 1+=recoil)
int _id;             // Unique identifier
int _tag;            // Material tag
Real _Ef;            // Final energy threshold
StateType _state;    // Current state (MOVING, REPLACEMENT, etc.)
```

**Ion States:**
- `MOVING` - Ion is still being tracked
- `REPLACEMENT` - Same element, energy below displacement threshold
- `SUBSTITUTIONAL` - Different element, energy below displacement threshold
- `INTERSTITIAL` - No recoil spawned, energy below final threshold
- `LOST` - Ion has left the sample
- `VACANCY` - Vacancy created

**Derived Classes:**
- `IonMDTag` - Includes MD energy gap tracking
- `IonClock` - Includes time tracking

#### MaterialBase (material.h:34)
Defines material composition and calculates stopping powers.

**Usage Pattern:**
```cpp
// Create material with density in g/cm³
MaterialBase * material = new MaterialBase(simconf, 8.92);

// Add elements
Element element;
element._Z = 29;        // Copper
element._m = 63.546;    // Atomic mass
element._t = 1.0;       // Stoichiometric fraction
material->_element.push_back(element);

// Must call prepare() after all elements are added
material->prepare();
```

**Key Members:**
```cpp
Real _rho;                    // Density in g/cm³
std::vector<Element> _element; // Elements in material
Real _am, _az;                // Average mass and atomic number
Real _arho;                   // Atomic density in atoms/Ang³
```

**Methods:**
- `prepare()` - Normalizes stoichiometry and computes averages (REQUIRED after adding elements)
- `average(const IonBase* pka)` - Computes PKA-dependent averages
- `getrstop(const IonBase* pka)` - Returns stopping power in eV/Angstrom

#### Element (element.h:29)
Represents a chemical element in a material.

**Key Members:**
```cpp
int _Z;          // Atomic number
Real _m;         // Atomic mass in amu
Real _t;         // Stoichiometric fraction (relative amount)
Real _Edisp;     // Displacement energy in eV (default: 25 eV)
Real _Elbind;    // Lattice binding energy in eV (default: 3 eV)
```

#### SampleBase (sample.h:33)
Abstract base class defining sample geometry.

**Key Members:**
```cpp
std::vector<MaterialBase*> material;  // Materials in the sample
Real w[3];                             // Simulation volume dimensions
sampleBoundary bc[3];                  // Boundary conditions per axis
```

**Boundary Conditions:**
- `PBC` - Periodic boundary conditions
- `INF` - Infinite (no boundaries)
- `CUT` - Cut off cascades at boundaries

**Virtual Methods to Override:**
```cpp
// Returns material at given position (MUST implement)
virtual MaterialBase* lookupMaterial(Point& pos) = 0;

// Returns distance to next material boundary (optional optimization)
virtual Real rangeMaterial(Point& pos, Point& dir);

// Called for averaging/statistics (optional)
virtual void averages(const IonBase* pka);
```

**Built-in Sample Geometries:**

1. **SampleLayers** (sample_layers.h:30) - Planar layered materials
   ```cpp
   SampleLayers* sample = new SampleLayers(thickness, width, height);
   sample->material.push_back(material1);
   sample->layerThickness.push_back(100.0);  // Angstroms
   sample->material.push_back(material2);
   sample->layerThickness.push_back(200.0);
   ```

2. **SampleSolid** (sample_solid.h) - Single uniform material
   ```cpp
   SampleSolid* sample = new SampleSolid(width, height, depth);
   sample->material.push_back(material);
   ```

3. **SampleWire** (sample_wire.h) - Cylindrical wire
4. **SampleBurriedWire** (sample_burried_wire.h) - Wire in matrix
5. **SampleClusters** (sample_clusters.h) - Spherical cluster inclusions
6. **SampleDynamic** (sample_dynamic.h) - Time-dependent materials

#### TrimBase (trim.h:36)
Base class for TRIM simulation runners. This is the main simulation engine.

**Constructor:**
```cpp
TrimBase(SimconfType* simconf, SampleBase* sample);
```

**Key Methods:**
```cpp
// Run simulation for a single PKA, generate recoils
void trim(IonBase* pka, std::queue<IonBase*>& recoils);

// Set output file base name
void setBaseName(const std::string& name);

// Write simulation results (override in subclasses)
virtual void writeOutput();
```

**Scattering Potentials:**
```cpp
enum Potential { UNIVERSAL, MOLIERE, CKR };
Potential _potential;  // Default: UNIVERSAL
```

**Virtual Hooks for Customization:**

All these methods are called at specific points during the cascade simulation:

```cpp
// Controls which recoils to follow (default: E > 12 eV)
virtual bool followRecoil();

// Called when a vacancy is created
virtual void vacancyCreation();

// Called when a replacement collision occurs
virtual void replacementCollision();

// Called to check PKA state
virtual void checkPKAState();

// Called when recoil energy is dissipated (phonons)
virtual void dissipateRecoilEnergy();
```

**Built-in TrimBase Subclasses:**

1. **TrimPrimaries** (trim.h:113) - Only follows primary knock-ons (gen 0)
2. **TrimRecoils** (trim.h:127) - Follows primaries and first-generation recoils
3. **TrimHistory** (trim.h:139) - Records full trajectory history
4. **TrimDefectLog** (trim.h:159) - Logs vacancy/interstitial creation to stream
5. **TrimVacMap** (trim.h:180) - Maps vacancy creation spatially
6. **TrimPhononOut** (trim.h:206) - Outputs phonon energy losses

**Example Usage:**
```cpp
// Initialize simulation
SimconfType* simconf = new SimconfType(seed);
simconf->tmin = 0.2;  // Minimum tracking energy

SampleSolid* sample = new SampleSolid(100.0, 100.0, 100.0);
sample->material.push_back(material);

TrimBase* trim = new TrimBase(simconf, sample);

// Create PKA
IonBase* pka = new IonBase();
pka->_Z = 29;           // Copper
pka->_m = 63.546;
pka->_E = 150000.0;     // 150 keV
pka->_gen = 0;
pka->_dir = Point(1.0, 0.0, 0.0);  // x direction
pka->_pos = Point(0.0, 50.0, 50.0);
pka->setEf();  // Set final energy threshold

// Run simulation
std::queue<IonBase*> recoils;
recoils.push(pka);

while (!recoils.empty()) {
    IonBase* ion = recoils.front();
    recoils.pop();
    sample->averages(ion);
    trim->trim(ion, recoils);
    delete ion;
}
```

---

## Using runmytrim

The `runmytrim` executable is the main command-line tool for running TRIM simulations using JSON input files.

### Basic Usage

```bash
./runmytrim < input.json
```

The executable reads a JSON input file from standard input and performs the simulation.

### JSON Input Format

#### Complete Example

```json
{
  "mytrim": {
    "options": {
      "seed": 12345,
      "threads": 4,
      "scale": 10
    },
    "ion": {
      "Z": 29,
      "mass": 63.546,
      "energy": 150000,
      "number": 1000
    },
    "sample": {
      "layers": [
        {
          "thickness": 1000,
          "rho": 8.92,
          "elements": [
            {
              "Z": 29,
              "mass": 63.546,
              "fraction": 1.0,
              "edisp": 25.0,
              "elbind": 3.0
            }
          ]
        }
      ]
    },
    "output": {
      "base": "simulation",
      "type": "vaccount",
      "primaries_only": false
    }
  }
}
```

### Input Sections

#### 1. Top-Level Block: `mytrim`

The entire input must be wrapped in a `mytrim` object:
```json
{
  "mytrim": {
    // All configuration goes here
  }
}
```

#### 2. Options Block (optional)

```json
"options": {
  "seed": <integer>,       // Random seed (optional, uses /dev/random if not specified)
  "threads": <integer>,    // Number of parallel threads (optional, default: 1)
  "scale": <number>        // Length scale (optional, default: 1 for Angstroms)
}
```

**Length Scale Values:**
- `1` = Angstroms (default)
- `10` = Nanometers
- `10000` = Micrometers

#### 3. Ion Block (required)

Defines the incident ion species and energy:

```json
"ion": {
  "Z": <integer>,          // Atomic number (REQUIRED)
  "mass": <number>,        // Atomic mass in amu (REQUIRED)
  "energy": <number>,      // Ion energy in eV (REQUIRED)
  "number": <integer>,     // Number of ions to simulate (REQUIRED)
  "final_energy": <number> // Final energy threshold in eV (optional)
}
```

**Example: 150 keV Copper ions**
```json
"ion": {
  "Z": 29,
  "mass": 63.546,
  "energy": 150000,
  "number": 10000
}
```

#### 4. Sample Block (required)

Defines the target material structure using layers:

```json
"sample": {
  "layers": [
    {
      "thickness": <number>,     // Layer thickness (in length scale units, REQUIRED)
      "rho": <number>,           // Density in g/cm³ (REQUIRED)
      "elements": [              // Array of elements (REQUIRED)
        {
          "Z": <integer>,        // Atomic number (REQUIRED)
          "mass": <number>,      // Atomic mass in amu (REQUIRED)
          "fraction": <number>,  // Stoichiometric fraction (REQUIRED)
          "edisp": <number>,     // Displacement energy in eV (optional, default: 25)
          "elbind": <number>     // Lattice binding energy in eV (optional, default: 3)
        }
      ]
    }
  ]
}
```

**Single Element Example (Cu):**
```json
"sample": {
  "layers": [
    {
      "thickness": 1000,
      "rho": 8.92,
      "elements": [
        { "Z": 29, "mass": 63.546, "fraction": 1.0 }
      ]
    }
  ]
}
```

**Compound Material Example (UO₂):**
```json
"sample": {
  "layers": [
    {
      "thickness": 10000,
      "rho": 10.97,
      "elements": [
        { "Z": 92, "mass": 238.03, "fraction": 1.0 },
        { "Z": 8, "mass": 15.999, "fraction": 2.0 }
      ]
    }
  ]
}
```

**Multi-Layer Example (Ag on Cu):**
```json
"sample": {
  "layers": [
    {
      "thickness": 100,
      "rho": 10.5,
      "elements": [
        { "Z": 47, "mass": 107.868, "fraction": 1.0 }
      ]
    },
    {
      "thickness": 1000,
      "rho": 8.92,
      "elements": [
        { "Z": 29, "mass": 63.546, "fraction": 1.0 }
      ]
    }
  ]
}
```

#### 5. Output Block (required)

Controls simulation output:

```json
"output": {
  "base": <string>,              // Output filename base (REQUIRED)
  "type": <string>,              // Output type (REQUIRED)
  "primaries_only": <boolean>    // Only track primary ions (optional, default: false)
}
```

**Output Types:**

1. **`vaccount`** - Vacancy counting
   - Produces: `<base>_vac.dat`
   - Format: `x vacancy_count replacement_count`
   - Each line represents counts per Angstrom depth

2. **`vacenergycount`** - Vacancy counting with energy binning
   - Produces: `<base>_vac.dat`
   - Bins vacancies by both depth and log₁₀(energy)

3. **`range`** - Ion range calculation
   - Produces: `<base>_ranges.dat`
   - Histogram of final ion positions
   - Separate column for each element (by Z)

**primaries_only Option:**
- `true` - Only track primary knock-ons (PKA), ignore recoils
- `false` - Track full cascade including all recoils (default)

### Complete Working Examples

#### Example 1: Self-Ion Irradiation
```json
{
  "mytrim": {
    "options": { "seed": 2344 },
    "ion": { "Z": 29, "mass": 63.546, "energy": 150000, "number": 1000 },
    "sample": {
      "layers": [
        { "thickness": 1000, "rho": 8.92, "elements": [
          { "Z": 29, "mass": 63.546, "fraction": 1 }
        ]}
      ]
    },
    "output": { "base": "cu_on_cu", "type": "vaccount" }
  }
}
```

#### Example 2: Heavy Ion on Nuclear Fuel
```json
{
  "mytrim": {
    "options": { "seed": 2344, "scale": 10, "threads": 8 },
    "ion": { "Z": 54, "mass": 131.904, "energy": 10000000, "number": 100 },
    "sample": {
      "layers": [
        { "thickness": 10000, "rho": 10.97, "elements": [
          { "Z": 92, "mass": 238.03, "fraction": 1 },
          { "Z": 8, "mass": 15.999, "fraction": 2 }
        ]}
      ]
    },
    "output": { "base": "xe_on_uo2", "type": "range" }
  }
}
```

### Understanding Output Files

#### Vacancy Count Output (`_vac.dat`)
```
# x(Angstrom)  vacancies  replacements
0 0 0
1 5 2
2 18 7
3 42 15
...
```

Each line shows:
- Depth in Angstroms
- Number of vacancies created at that depth
- Number of replacement collisions at that depth

#### Range Output (`_ranges.dat`)
```
#x Z29 Z47
0 0 0
10 2 0
20 8 1
30 15 3
...
```

Each line shows:
- x position (depth)
- Count histogram for each element that stopped in the sample

### Multi-Threading

The `runmytrim` executable supports multi-threaded execution:

```json
"options": { "threads": 8 }
```

**Threading Architecture:**
- Each thread has its own `SimconfType`, `SampleBase`, and `TrimBase` instance
- PKAs are distributed round-robin across threads
- Results are joined at the end
- Scales well up to the number of physical cores

**Thread Safety Notes:**
- Each thread uses a different RNG seed (derived from master seed)
- No synchronization needed during simulation
- Joining happens only at the end

---

## Using runstopping

The `runstopping` executable calculates electronic stopping power as a function of ion energy.

### Basic Usage

```bash
./runstopping < input.json
```

### JSON Input Format

#### Complete Example

```json
{
  "stopping": {
    "ion": {
      "Z": 1,
      "mass": 1.008,
      "energy": { "begin": 10000, "end": 10000000, "mult": 1.1 }
    },
    "material": {
      "rho": 7.8658,
      "elements": [
        { "Z": 26, "mass": 55.847, "fraction": 1.0 }
      ]
    }
  }
}
```

### Input Sections

#### 1. Top-Level Block: `stopping`

```json
{
  "stopping": {
    // All configuration goes here
  }
}
```

#### 2. Ion Block (required)

```json
"ion": {
  "Z": <integer>,      // Atomic number (REQUIRED)
  "mass": <number>,    // Atomic mass in amu (REQUIRED)
  "energy": <...>      // Energy specification (REQUIRED)
}
```

**Energy Specification - Three Options:**

**Option 1: Single Energy Value**
```json
"energy": 1000000
```

**Option 2: Array of Specific Energies**
```json
"energy": [10, 100, 1000, 10000, 100000, 1000000]
```

**Option 3: Energy Range with Stepping**

Linear stepping:
```json
"energy": {
  "begin": 1000,
  "end": 100000,
  "step": 1000
}
```

Logarithmic stepping (recommended for wide ranges):
```json
"energy": {
  "begin": 10000,
  "end": 10000000,
  "mult": 1.1
}
```

For logarithmic stepping:
- `mult` must be > 1.0
- Each energy is: E(n+1) = E(n) × mult
- Provides even spacing on log scale

#### 3. Material Block (required)

```json
"material": {
  "rho": <number>,           // Density in g/cm³ (REQUIRED)
  "elements": [              // Array of elements (REQUIRED)
    {
      "Z": <integer>,        // Atomic number (REQUIRED)
      "mass": <number>,      // Atomic mass in amu (REQUIRED)
      "fraction": <number>   // Stoichiometric fraction (REQUIRED)
    }
  ]
}
```

### Output Format

The output is written to stdout in two-column format:

```
<energy_eV> <stopping_power_eV_per_Angstrom>
```

Example output:
```
10000 156.234
11000 165.891
12100 175.123
...
```

### Complete Working Examples

#### Example 1: Hydrogen in Iron
```json
{
  "stopping": {
    "ion": {
      "Z": 1,
      "mass": 1.008,
      "energy": { "begin": 1000, "end": 1000000, "mult": 1.2 }
    },
    "material": {
      "rho": 7.8658,
      "elements": [
        { "Z": 26, "mass": 55.847, "fraction": 1.0 }
      ]
    }
  }
}
```

#### Example 2: Helium in Silicon
```json
{
  "stopping": {
    "ion": {
      "Z": 2,
      "mass": 4.003,
      "energy": [100, 500, 1000, 5000, 10000, 50000, 100000]
    },
    "material": {
      "rho": 2.329,
      "elements": [
        { "Z": 14, "mass": 28.085, "fraction": 1.0 }
      ]
    }
  }
}
```

#### Example 3: Copper in Tungsten
```json
{
  "stopping": {
    "ion": {
      "Z": 29,
      "mass": 63.546,
      "energy": { "begin": 10000, "end": 10000000, "step": 10000 }
    },
    "material": {
      "rho": 19.25,
      "elements": [
        { "Z": 74, "mass": 183.84, "fraction": 1.0 }
      ]
    }
  }
}
```

### Usage Tips

1. **Use logarithmic stepping for wide energy ranges** - Provides better resolution where stopping power changes rapidly
2. **Output can be piped** - Redirect to file or pipe to plotting tools:
   ```bash
   ./runstopping < input.json > stopping_data.txt
   ./runstopping < input.json | gnuplot -e "plot '-' with lines"
   ```
3. **Units are always eV and Angstroms** - No length scale option for runstopping

---

## Customizing MyTRIM

MyTRIM is designed to be extended through subclassing. The most common customization points are `TrimBase` and `SampleBase`.

### Creating Custom TRIM Simulations

Subclass `TrimBase` or `ThreadedTrimBase` to customize simulation behavior.

#### Example 1: Custom Vacancy Counter

```cpp
#include "trim.h"
#include <fstream>
#include <vector>

class TrimVacCount : public ThreadedTrimBase
{
public:
  TrimVacCount(SimconfType * simconf, SampleBase * sample)
    : ThreadedTrimBase(simconf, sample)
  {
  }

protected:
  // Called when a vacancy is created
  virtual void vacancyCreation() override
  {
    _simconf->vacancies_created++;

    // Bin by x position
    int x = _recoil->_pos(0);
    if (x >= 0) {
      if (x >= int(_vac_bin.size()))
        _vac_bin.resize(x + 1, 0);
      _vac_bin[x]++;
    }
  }

  // Called when replacement collision occurs
  virtual void replacementCollision() override
  {
    int x = _recoil->_pos(0);
    if (x >= 0) {
      if (x >= int(_repl_bin.size()))
        _repl_bin.resize(x + 1, 0);
      _repl_bin[x]++;
    }
  }

  // Join results from another thread
  virtual void threadJoin(const ThreadedTrimBase & ttb) override
  {
    const TrimVacCount & tvc = static_cast<const TrimVacCount &>(ttb);

    // Merge histogram data
    _vac_bin.resize(std::max(_vac_bin.size(), tvc._vac_bin.size()));
    for (unsigned int x = 0; x < tvc._vac_bin.size(); ++x)
      _vac_bin[x] += tvc._vac_bin[x];

    _repl_bin.resize(std::max(_repl_bin.size(), tvc._repl_bin.size()));
    for (unsigned int x = 0; x < tvc._repl_bin.size(); ++x)
      _repl_bin[x] += tvc._repl_bin[x];
  }

  // Write results to file
  virtual void writeOutput() override
  {
    std::ofstream out((_base_name + "_vac.dat").c_str());

    const unsigned int size = std::max(_vac_bin.size(), _repl_bin.size());
    _vac_bin.resize(size);
    _repl_bin.resize(size);

    for (unsigned int x = 0; x < size; ++x)
      out << x << ' ' << _vac_bin[x] << ' ' << _repl_bin[x] << '\n';
  }

private:
  std::vector<unsigned int> _vac_bin;
  std::vector<unsigned int> _repl_bin;
};
```

#### Example 2: Ion Range Tracking

```cpp
class TrimRange : public ThreadedTrimBase
{
public:
  TrimRange(SimconfType * simconf, SampleBase * sample)
    : ThreadedTrimBase(simconf, sample), _range(112)  // Space for all Z
  {
  }

protected:
  // Only follow primary ions, not recoils
  virtual bool followRecoil() override
  {
    return (_recoil->_gen < 1);
  }

  // Record final position when ion stops
  virtual void dissipateRecoilEnergy() override
  {
    _range[_recoil->_Z].push_back(_recoil->_pos(0));
  }

  // Vacancy calculation (Kinchin-Pease model)
  virtual void vacancyCreation() override
  {
    const Real Ed = _element->_Edisp;
    const Real ed = 0.0115 * std::pow(_recoil->_Z, -7.0/3.0) * _recoil->_E;
    const Real kd = 0.1337 * std::pow(_recoil->_Z, 2.0/3.0) / std::sqrt(_recoil->_m);
    const Real g = 3.4008 * std::pow(ed, 1.0/6.0) + 0.40244 * std::pow(ed, 3.0/4.0) + ed;
    const Real Ev = _recoil->_E / (1.0 + kd * g);

    if (Ev < Ed) return;

    if (Ev >= Ed / 0.4)
      _simconf->vacancies_created += Ev * 0.4 / Ed;
    else
      _simconf->vacancies_created++;
  }

  virtual void threadJoin(const ThreadedTrimBase & ttb) override
  {
    const TrimRange & tr = static_cast<const TrimRange &>(ttb);

    for (unsigned int Z = 0; Z < _range.size(); ++Z)
      _range[Z].insert(_range[Z].end(), tr._range[Z].begin(), tr._range[Z].end());
  }

  virtual void writeOutput() override
  {
    // Build histograms and write to file
    // (implementation details omitted for brevity)
  }

private:
  std::vector<std::vector<Real>> _range;  // Range data per element
};
```

#### Key Virtual Methods to Override

**`bool followRecoil()`**
- Controls which recoils are tracked in the cascade
- Return `true` to continue tracking, `false` to stop
- Access `_recoil` to inspect the recoil properties
- Default: follows if E > 12 eV

```cpp
// Only track primaries
virtual bool followRecoil() override {
  return (_recoil->_gen == 0);
}

// Track up to 2nd generation
virtual bool followRecoil() override {
  return (_recoil->_gen < 2);
}

// Track only specific element
virtual bool followRecoil() override {
  return (_recoil->_Z == 54);  // Only Xe
}
```

**`void vacancyCreation()`**
- Called when an atom is displaced from its lattice site
- Access `_recoil` for the displaced atom
- Access `_element` for the target element
- Use to count/record/bin vacancy events

**`void replacementCollision()`**
- Called when same-species atom replaces another below Edisp
- Typically used for tracking replacement sequences

**`void checkPKAState()`**
- Called when PKA state changes
- Access `_pka` for the ion
- Useful for recording final ion states/positions

**`void dissipateRecoilEnergy()`**
- Called when recoil energy is dissipated as phonons
- Ion comes to rest without creating vacancy
- Access `_recoil` for final position/energy

**`void writeOutput()`**
- Called at end of simulation to write results
- Use `_base_name` for output filenames

#### ThreadedTrimBase

For multi-threaded simulations, use `ThreadedTrimBase`:

```cpp
class MyCustomTrim : public ThreadedTrimBase
{
public:
  MyCustomTrim(SimconfType * simconf, SampleBase * sample)
    : ThreadedTrimBase(simconf, sample)
  {
  }

protected:
  // Must implement threadJoin for multi-threading
  virtual void threadJoin(const ThreadedTrimBase & ttb) override
  {
    const MyCustomTrim & other = static_cast<const MyCustomTrim &>(ttb);
    // Merge data from other thread into this thread
  }

  // Other virtual methods as needed...
};
```

The `_primaries_only` member can be set from JSON input:
```cpp
bool _primaries_only;  // Controlled by JSON "primaries_only" option
```

### Creating Custom Sample Geometries

Subclass `SampleBase` to define custom material geometries.

#### Example: Custom Spherical Sample

```cpp
#include "sample.h"

class SampleSphere : public SampleBase
{
public:
  SampleSphere(Real radius, Real box_size)
    : SampleBase(box_size, box_size, box_size), _radius(radius)
  {
    bc[0] = bc[1] = bc[2] = CUT;  // Cut off at boundaries
  }

  // Must implement: return material at position
  virtual MaterialBase * lookupMaterial(Point & pos) override
  {
    // Calculate distance from center
    Point center(w[0]/2.0, w[1]/2.0, w[2]/2.0);
    Real dist = (pos - center).norm();

    if (dist < _radius)
      return material[0];  // Inside sphere
    else
      return nullptr;      // Outside (vacuum)
  }

  // Optional: optimize by calculating intersection distance
  virtual Real rangeMaterial(Point & pos, Point & dir) override
  {
    Point center(w[0]/2.0, w[1]/2.0, w[2]/2.0);
    Point rel_pos = pos - center;

    // Ray-sphere intersection
    Real a = dir.norm_sq();
    Real b = 2.0 * (rel_pos * dir);
    Real c = rel_pos.norm_sq() - _radius * _radius;
    Real discriminant = b*b - 4*a*c;

    if (discriminant < 0)
      return 1e10;  // No intersection

    Real t1 = (-b - std::sqrt(discriminant)) / (2*a);
    Real t2 = (-b + std::sqrt(discriminant)) / (2*a);

    if (t1 > 0) return t1;
    if (t2 > 0) return t2;
    return 1e10;
  }

private:
  Real _radius;
};
```

#### Example: Multi-Layer with Varying Composition

```cpp
class SampleGradient : public SampleBase
{
public:
  SampleGradient(Real thickness, Real width, Real height, int n_layers)
    : SampleBase(thickness, width, height), _n_layers(n_layers)
  {
    bc[0] = CUT;  // Cut in depth direction
    bc[1] = bc[2] = PBC;  // Periodic in lateral directions

    // Create materials with varying composition
    for (int i = 0; i < n_layers; ++i) {
      Real fraction = Real(i) / (n_layers - 1);
      MaterialBase * mat = createGradientMaterial(fraction);
      material.push_back(mat);
    }
  }

  virtual MaterialBase * lookupMaterial(Point & pos) override
  {
    if (pos(0) < 0 || pos(0) > w[0])
      return nullptr;

    int layer = int(pos(0) * _n_layers / w[0]);
    if (layer >= _n_layers) layer = _n_layers - 1;

    return material[layer];
  }

  virtual Real rangeMaterial(Point & pos, Point & dir) override
  {
    // Distance to next layer boundary
    Real layer_thickness = w[0] / _n_layers;
    int current_layer = int(pos(0) / layer_thickness);
    Real next_boundary = (current_layer + 1) * layer_thickness;

    if (dir(0) > 0)
      return (next_boundary - pos(0)) / dir(0);
    else if (dir(0) < 0)
      return -(pos(0) - current_layer * layer_thickness) / dir(0);
    else
      return 1e10;  // Parallel to layers
  }

private:
  int _n_layers;

  MaterialBase * createGradientMaterial(Real fraction)
  {
    // Create material with composition based on fraction
    // e.g., Cu(1-f) + Ni(f)
    MaterialBase * mat = new MaterialBase(_simconf, 8.9);

    Element cu, ni;
    cu._Z = 29; cu._m = 63.546; cu._t = 1.0 - fraction;
    ni._Z = 28; ni._m = 58.693; ni._t = fraction;

    mat->_element.push_back(cu);
    mat->_element.push_back(ni);
    mat->prepare();

    return mat;
  }
};
```

#### Key Points for Custom Geometries

1. **`lookupMaterial()` must be fast** - Called many times per ion
2. **Return `nullptr` for vacuum** - Ion will be marked as LOST
3. **`rangeMaterial()` is optional but recommended** - Speeds up simulation by reducing material lookups
4. **Set boundary conditions appropriately**:
   - `PBC` for infinite periodic structures
   - `INF` for effectively infinite samples
   - `CUT` to terminate ions at boundaries
5. **Materials are not thread-safe** - Each thread should have its own material instances

### Using Custom Classes with runmytrim

To use custom TRIM classes with `runmytrim`:

1. Add your class header to `apps/include/`
2. Add implementation to `apps/src/`
3. Modify `runmytrim.C` to recognize your output type:

```cpp
if (output_type == "myoutput")
  for (auto & td : thread_data)
    td._trim = new MyCustomTrim(&(td._simconf), td._sample);
```

4. Rebuild:
```bash
cd build
make
```

---

## Building and Installation

### Prerequisites

- C++11 compatible compiler (GCC 4.8+, Clang 3.3+)
- CMake 2.8+
- JsonCpp library (for `runmytrim` and `runstopping` only)

### Installing JsonCpp

Required for the JSON-based executables:

```bash
git clone https://github.com/open-source-parsers/jsoncpp.git
cd jsoncpp
mkdir build && cd build
cmake -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF -G "Unix Makefiles" ..
sudo make install
```

### Building MyTRIM

```bash
# Clone repository
git clone https://github.com/idaholab/mytrim.git
cd mytrim

# Standard build
mkdir build && cd build
cmake ..
make

# Install (optional)
sudo make install
```

### Build Options

**Debug build:**
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```
- Enables `-g -O1 -Wall -Wextra -Werror -pedantic`
- Useful for development

**Performance profiling with gperftools:**
```bash
cmake -DENABLE_GPERFTOOLS=ON ..
make
```

**LLVM profiling instrumentation:**
```bash
cmake -DENABLE_LLVM_PROFILE=ON ..
make
```

### Installation Paths

Default installation paths (with `-DCMAKE_INSTALL_PREFIX=/usr/local`):

- Library: `/usr/local/lib/libmytrim.a`
- Headers: `/usr/local/include/mytrim/`
- Data files: `/usr/local/share/mytrim/`

### Testing

```bash
cd tests
./run_tests.sh
```

Individual tests:
```bash
cd tests/json
./test.sh
```

### Data Files

MyTRIM requires data files from SRIM in the `data/` directory:
- `SCOEF.95A` - Stopping coefficients (part A)
- `SCOEF.95B` - Stopping coefficients (part B)
- `SLFCTR.dat` - Screening length factors
- `SNUC03.dat` - Nuclear stopping data
- `ELNAME.dat` - Element names

These are automatically installed to the data directory during `make install`.

### Environment Variables

If running from build directory without installing:
```bash
export MYTRIM_DATADIR=/path/to/mytrim/data
```

### Integration with MOOSE

To build as part of MOOSE framework:

1. The `MYTRIM_ENABLED` preprocessor flag is set by MOOSE
2. Uses MOOSE types: `Real`, `Point`, `MooseError`
3. Integrated via the Magpie application

---

## Additional Resources

- **Repository:** https://github.com/idaholab/mytrim
- **License:** LGPL 2.1
- **SRIM Website:** http://www.srim.org
- **Wikipedia:** https://en.wikipedia.org/wiki/Stopping_and_Range_of_Ions_in_Matter

## Citation

If you use MyTRIM in your research, please cite:

*MyTRIM - A three dimensional binary collision Monte Carlo library*
Daniel Schwen, Idaho National Laboratory

---

## Troubleshooting

### Common Issues

**1. JsonCpp not found**
```
CMake Error: Could not find JsonCpp
```
Install JsonCpp (see prerequisites) or the executables won't be built.

**2. Data files not found**
```
ERROR: Unable to open data file
```
Set `MYTRIM_DATADIR` environment variable or ensure files are in `/usr/local/share/mytrim/`.

**3. Segmentation fault in simulation**
```
Segmentation fault (core dumped)
```
Most common causes:
- Forgot to call `material->prepare()` after adding elements
- `lookupMaterial()` returning invalid pointer
- Accessing `_recoil` or `_pka` when it's null

**4. Test failures**
```
Difference in out.xxx
```
May be due to:
- Different random number generation (check MYTRIM_SEED)
- Different compiler optimizations
- Numerical precision differences

**5. Multi-threading issues**
- Ensure each thread has its own `SimconfType` instance with different seed
- Materials should not be shared between threads during simulation
- Only join data after all threads complete

### Performance Tips

1. **Use appropriate length scales** - Reduces numerical errors
2. **Enable optimizations** - Use Release build (`-O3`)
3. **Adjust `tmin`** - Higher values = faster but less accurate
4. **Use `primaries_only`** - Much faster if recoils not needed
5. **Optimize `rangeMaterial()`** - Reduces material lookups
6. **Use multi-threading** - Scales linearly with cores
