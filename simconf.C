#include "simconf.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <cstdlib>

namespace MyTRIM_NS {
  SimconfType *simconf;
}

using namespace MyTRIM_NS;

SimconfType::SimconfType(Real _alfa) :
    scoef(_rows),
    _data_dir(std::getenv("MYTRIM_DATADIR") ? std::getenv("MYTRIM_DATADIR") : MYTRIM_DATA_DIR)
{
  ed = 25.0; // displacement energy
  alfa = _alfa; // angle of incidence (degrees)
  alpha = alfa * M_PI / 180.0;
  tmin = 1.0; //max impact parameter set by min. transferred energy
  //tmin = 5.0; //max impact parameter set by min. transferred energy
  tau = 0.0;
  da = 3.0; // angular grid for transmitted ions
  cw = 0.001; // channel width 1% of layer thickness

  // output full trajectories
  fullTraj = false;

  // set global ion id to zero (will be incremented for each new projectile)
  id = 0;

  // initialize global statistics
  vacancies_created = 0;
  EelTotal = 0.0;
  EnucTotal = 0.0;

  // read data tables
  readDataFiles();
}

void
SimconfType::readDataFiles()
{
  std::string fname;
  int dummy;

  fname = _data_dir + "/SNUC03.dat";
  std::ifstream snuc03(fname.c_str());
  if (!snuc03)
    fileReadError(fname);

  for (int i = 0; i < 92; ++i)
    for (int j = i; j < 92; j++)
    {
      snuc03 >> dummy >> dummy
             >> snuc[j][i][0] >> snuc[j][i][1] >> snuc[j][i][2] >> snuc[j][i][3];

      if (!snuc03)
        fileReadError("contents of SNUC03.dat");

      for (int n = 0; n < 4; n++)
        snuc[i][j][n] = snuc[j][i][n];
    }


  fname = _data_dir + "/SCOEF.95A";
  std::ifstream scoef95a(fname.c_str());
  if (!scoef95a)
    fileReadError(fname);

  skipLine(scoef95a); // header
  skipLine(scoef95a); // header

  fname = _data_dir + "/SCOEF.95B";
  std::ifstream scoef95b(fname.c_str());
  if (!scoef95b)
    fileReadError(fname);
  skipLine(scoef95b); // header
  skipLine(scoef95b); // header


  fname = _data_dir + "/SLFCTR.dat";
  std::ifstream slfctr(fname.c_str());
  if (!slfctr)
    fileReadError(fname);
  skipLine(slfctr); // header

  fname = _data_dir + "/ELNAME.dat";
  std::ifstream elname(fname.c_str());
  if (!elname)
    fileReadError(fname);

  for (unsigned int i = 0; i < _rows; ++i)
  {
    scoef[i].read95A(scoef95a);
    scoef[i].read95B(scoef95b);
    scoef[i].readSlfctr(slfctr);
    scoef[i].readElname(elname);
  }

  // read the last line from the 95B data which contains energy range data
  scoeflast.read95B(scoef95b);
}

void
SimconfType::fileReadError(const std::string & path)
{
#ifdef MYTRIM_ENABLED
    mooseError("Error reading " << path);
#else
    std::cerr << "Error reading " << path << std::endl;
    exit(1);
#endif
}

void
SimconfType::skipLine(std::ifstream & sf)
{
  std::string s;
  if (!std::getline(sf, s))
  {
    #ifdef MYTRIM_ENABLED
      mooseError("Error reading file");
    #else
      std::cerr << "Error reading file" << std::endl;
      exit(1);
    #endif
  }
}


SimconfType::ScoefLine::ScoefLine() :
    mm1(0.0),
    m1(0.0),
    mnat(0.0),
    rho(0.0),
    atrho(0.0),
    vfermi(0.0),
    heat(0.0),
    lfctr(0.0),
    pcoef(8),
    ehigh(4),
    screen(19),
    fermicorr(15)
{
}

void
SimconfType::ScoefLine::read95A(std::ifstream & scoef95a)
{
  int dummy;
  scoef95a >> dummy
           >> mm1 >> m1 >> mnat
           >> rho >> atrho >> vfermi >> heat;

  if (!scoef95a)
    SimconfType::fileReadError("contents of SCOEF.95A");

  for (unsigned int j = 0; j < 8; ++j)
    if (!(scoef95a >> pcoef[j]))
      SimconfType::fileReadError("proton coefficient in SCOEF.95B");
}

void
SimconfType::ScoefLine::read95B(std::ifstream & scoef95b)
{
  for (unsigned int j = 0; j < 4; ++j)
    if (!(scoef95b >> ehigh[j]))
      SimconfType::fileReadError("high energy coefficient in SCOEF.95B");
  for (unsigned int j = 0; j < 19; ++j)
    if (!(scoef95b >> screen[j]))
      SimconfType::fileReadError("screening data in SCOEF.95B");
  for (unsigned int j = 0; j < 15; ++j)
    if (!(scoef95b >> fermicorr[j]))
      SimconfType::fileReadError("fermi correction in SCOEF.95B");
}

void
SimconfType::ScoefLine::readSlfctr(std::ifstream & slfctr)
{
  int dummy;
  if (!(slfctr >> dummy >> lfctr))
    SimconfType::fileReadError("contents of SLFCTR.dat");
}

void
SimconfType::ScoefLine::readElname(std::ifstream & elname)
{
  int dummy;
  if (!(elname >> dummy >> sym >> name))
    SimconfType::fileReadError("contents of ELNAME.dat");
}
