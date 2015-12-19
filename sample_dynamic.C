#include "sample_dynamic.h"
#include "simconf.h"
#include "element.h"
#include <iostream>
#include <cmath>

using namespace MyTRIM_NS;

sampleDynamic::sampleDynamic(simconfType * simconf_, Real x, Real y, Real z):
    sampleLayers(x, y, z),
    simconf(simconf_)
{
  bc[0] = CUT;
  bc[1] = PBC;
  bc[2] = PBC;
}

void sampleDynamic::averages(const ionBase *_pka)
{
  // remember pka, we do not calculate averages right now, but on demand!
  pka = _pka;

  // reset update status, to trigger on-demand updates
  int i = material.size();
  while (i > 0) material[--i]->dirty = true;
}

materialBase*  sampleDynamic::lookupMaterial(Real* pos)
{
  //std::cout << "lookuplayer" << std::endl;
  materialBase* m = material[lookupLayer(pos)];
  //std::cout << m << ' ' << m->arho << ' ' << m->am << ' ' << m->az << ' ' << m->mu << std::endl;

  // on-demand update
  if (m->dirty)
  {
    //std::cout << "on demand aver" << std::endl;
    m->average(pka);
    //std::cout << m << ' ' << m->arho << ' ' << m->am << ' ' << m->az << ' ' << m->mu << std::endl;
  }
  return m;
}

void sampleDynamic::addAtomsToLayer(int layer, int n, int Z)
{
  unsigned int i, ne = material[layer]->element.size();
  Real rnorm = 0.0; // volume of one atom with the fractional composition of the layer at natural density

  // look which element in the layer we are modifying and take weighted average of 1/rho of elements
  for (i = 0; i < material[layer]->element.size(); ++i)
  {
    if (material[layer]->element[i]->z == Z) ne = i;
    rnorm += material[layer]->element[i]->t / simconf->scoef[material[layer]->element[i]->z-1].atrho;
  }

  // element not yet contained in layer
  if (ne == material[layer]->element.size())
  {
    elementBase* element = new elementBase;
    element->z = Z;
    element->m = simconf->scoef[Z-1].m1;
    element->t = 0.0;
    material[layer]->element.push_back(element);
  }

  // mass of layer (g/cm^3 * Ang^3 = 10^-24 g)
  //Real ml = material[layer]->rho * w[1] * w[2] * layerThickness[layer];

  // number of atoms in layer
  int nl    = material[layer]->arho * w[1] * w[2] * layerThickness[layer];

  // mass change of layer (g/mole = g/(0.6022*10^24))
  Real ma = 0.6022 * material[layer]->element[ne]->m * Real(n);

  // keep density consistent...
  layerThickness[layer] += ma / (material[layer]->rho * w[1] * w[2]);

  // change stoichiometry (arho 1/ang^3 * ang^3)
  if (material[layer]->element[ne]->t*nl + Real(n) < 0.0)
  {
    std::cout << "Crap, t*nl=" << material[layer]->element[ne]->t*nl << ", but n=" << n << std::endl;
    material[layer]->element[ne]->t = 0.0;
  }
  else
    material[layer]->element[ne]->t += Real(n)/nl; // sum t will be scaled to 1.0 in material->prepare()

  // mark as dirty, just in case, and re-prepare to update averages
  material[layer]->dirty = true;
  material[layer]->prepare();
}
