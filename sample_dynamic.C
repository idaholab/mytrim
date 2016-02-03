#include "sample_dynamic.h"
#include "simconf.h"
#include "element.h"
#include <iostream>
#include <cmath>

using namespace MyTRIM_NS;

SampleDynamic::SampleDynamic(SimconfType * simconf, Real x, Real y, Real z):
    SampleLayers(x, y, z),
    _simconf(simconf)
{
  bc[0] = CUT;
  bc[1] = PBC;
  bc[2] = PBC;
}

void
SampleDynamic::averages(const IonBase *_pka)
{
  // remember pka, we do not calculate averages right now, but on demand!
  pka = _pka;

  // reset update status, to trigger on-demand updates
  int i = material.size();
  while (i > 0) material[--i]->_dirty = true;
}

MaterialBase*
SampleDynamic::lookupMaterial(Point & pos)
{
  //std::cout << "lookuplayer" << std::endl;
  MaterialBase* m = material[lookupLayer(pos)];
  //std::cout << m << ' ' << m->_arho << ' ' << m->am << ' ' << m->az << ' ' << m->mu << std::endl;

  // on-demand update
  if (m->_dirty)
  {
    //std::cout << "on demand aver" << std::endl;
    m->average(pka);
    //std::cout << m << ' ' << m->_arho << ' ' << m->am << ' ' << m->az << ' ' << m->mu << std::endl;
  }
  return m;
}

void
SampleDynamic::addAtomsToLayer(int layer, int n, int Z)
{
  unsigned int i, ne = material[layer]->_element.size();
  Real rnorm = 0.0; // volume of one atom with the fractional composition of the layer at natural density

  // look which element in the layer we are modifying and take weighted average of 1/rho of elements
  for (i = 0; i < material[layer]->_element.size(); ++i)
  {
    if (material[layer]->_element[i]->_Z == Z) ne = i;
    rnorm += material[layer]->_element[i]->_t / _simconf->scoef[material[layer]->_element[i]->_Z-1].atrho;
  }

  // element not yet contained in layer
  if (ne == material[layer]->_element.size())
  {
    ElementBase* element = new ElementBase;
    element->_Z = Z;
    element->_m = _simconf->scoef[Z-1].m1;
    element->_t = 0.0;
    material[layer]->_element.push_back(element);
  }

  // mass of layer (g/cm^3 * Ang^3 = 10^-24 g)
  //Real ml = material[layer]->_rho * w[1] * w[2] * layerThickness[layer];

  // number of atoms in layer
  int nl    = material[layer]->_arho * w[1] * w[2] * layerThickness[layer];

  // mass change of layer (g/mole = g/(0.6022*10^24))
  Real ma = 0.6022 * material[layer]->_element[ne]->_m * Real(n);

  // keep density consistent...
  layerThickness[layer] += ma / (material[layer]->_rho * w[1] * w[2]);

  // change stoichiometry (_arho 1/ang^3 * ang^3)
  if (material[layer]->_element[ne]->_t*nl + Real(n) < 0.0)
  {
    std::cout << "Crap, t*nl=" << material[layer]->_element[ne]->_t*nl << ", but n=" << n << std::endl;
    material[layer]->_element[ne]->_t = 0.0;
  }
  else
    material[layer]->_element[ne]->_t += Real(n)/nl; // sum t will be scaled to 1.0 in material->prepare()

  // mark as _dirty, just in case, and re-prepare to update averages
  material[layer]->_dirty = true;
  material[layer]->prepare();
}
