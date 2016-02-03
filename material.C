#include "material.h"
#include "simconf.h"

#include "functions.h"

#include <cmath>
#include <iostream>

using namespace MyTRIM_NS;

MaterialBase::MaterialBase(SimconfType * simconf, Real rho) :
    _rho(rho),
    tag(-1),
    _dirty(true),
    _simconf(simconf)
{
}

void
MaterialBase::prepare()
{
  Real tt = 0.0;

  // get total stoichiometry
  for (unsigned int i = 0; i < _element.size(); ++i)
  {
    if (_element[i]->_t < 0.0)
      _element[i]->_t = 0.0;
    tt += _element[i]->_t;
  }

  // normalize relative probabilities to 1
  for (unsigned int i = 0; i < _element.size(); ++i)
    _element[i]->_t /= tt;

  // average
  _am = 0.0;
  _az = 0.0;
  for (unsigned int i = 0; i < _element.size(); ++i)
  {
    _am += _element[i]->_m * _element[i]->_t;
    _az += Real(_element[i]->_Z) * _element[i]->_t;
  }

  _arho = _rho * 0.6022 / _am; //[TRI00310] atoms/Ang^3
}

// make sure layers are prepare'd first!
void
MaterialBase::average(const IonBase * pka)
{
  mu = pka->_m / _am;

  // universal or firsov screening length
  a = .5292 * .8853 / (std::pow(Real(pka->_Z), 0.23) + std::pow(_az, 0.23));
  //a = .5292 * .8853 / std::pow(pow(Real(pka._Z), 0.5) + std::pow(_az, 0.5), 2.0/3.0);

  // mean flight path0
  f = a * _am / (_az * Real(pka->_Z) * 14.4 * (pka->_m + _am));
  //eps0 = e0 * f;
  epsdg = _simconf->tmin * f * std::pow(1.0 + mu, 2.0) / (4.0 * mu);

  // fd and kd determine how much recoil energy goes into el. loss and vaccancies
  fd = std::pow(0.01 * _az, -7.0 / 3.0);
  kd = std::pow(0.1334 * _az, 2.0 / 3.0) / std::sqrt(_am);

  for (unsigned int i = 0; i < _element.size(); ++i)
  {
    _element[i]->my = pka->_m / _element[i]->_m;
    _element[i]->ec = 4.0 * _element[i]->my / std::pow(1.0 + _element[i]->my, 2.0);
    _element[i]->ai = .5292 * .8853 / (std::pow(Real(pka->_Z), 0.23) + std::pow(_element[i]->_Z, 0.23));
    //ai = .5292 * .8853 / std::pow(pow(Real(pka._Z), 0.5) + std::pow(_element[i].z, 0.5), 2.0/3.0);
    _element[i]->fi = _element[i]->ai * _element[i]->_m /
                     (Real(pka->_Z) * Real(_element[i]->_Z) * 14.4 * (pka->_m + _element[i]->_m));
  }

  _dirty = false;
}

// make sure layers are prepare'd and averaged first!
Real
MaterialBase::getrstop(const IonBase * pka)
{
  Real se = 0.0;
  for (unsigned int i = 0; i < _element.size(); ++i)
    se += rstop(pka, _element[i]->_Z) * _element[i]->_t * _arho;

  return se;
}

Real
MaterialBase::rpstop(int z2p, Real e)
{
  Real pe, pe0, sl, sh, sp, velpwr;
  int z2 = z2p-1;
  // velocity proportional stopping below pe0
  pe0 = 25.0;
  pe = std::max(pe0, e);

  // pcoef indices are one less than in the fortran version!
  sl = (_simconf->pcoef[z2][0] * std::pow(pe, _simconf->pcoef[z2][1])) +
       (_simconf->pcoef[z2][2] * std::pow(pe, _simconf->pcoef[z2][3]));
  sh = _simconf->pcoef[z2][4] / std::pow(pe, _simconf->pcoef[z2][5]) *
       std::log(_simconf->pcoef[z2][6] / pe + _simconf->pcoef[z2][7] * pe);
  sp = sl * sh / (sl + sh);
  if (e <= pe0)
  {
    // velpwr is the power of velocity stopping below pe0
    if (z2p <= 6)
      velpwr = 0.25;
    else
      velpwr = 0.45;
    sp *= std::pow(e/pe0, velpwr);
  }
  return sp;
}

Real
MaterialBase::rstop(const IonBase * ion, int z2)
{
  Real e, vrmin, yrmin, v, vr, yr, vmin, m1;
  Real a, b, q, q1, l, l0, l1;
  Real zeta;
  int z1 = ion->_Z;
  Real fz1 = Real(z1), fz2 = Real(z2);
  Real eee, sp, power;
  Real se;

  // scoeff
  Real lfctr = _simconf->scoef[z1-1].lfctr;
  Real mm1 = _simconf->scoef[z1-1].mm1;
  Real vfermi = _simconf->scoef[z2-1].vfermi;
  //Real atrho = _simconf->scoef[z2-1].atrho;

  if (ion->_m == 0.0)
    m1 = mm1;
  else
    m1 = ion->_m;

  e = 0.001 * ion->_E / m1;

  if (z1 == 1)
  {
#ifdef MYTRIM_ENABLED
    mooseError("proton stopping not yet implemented!");
#else
    std::cerr << "proton stopping not yet implemented!\n";
    exit(1);
#endif
  }
  else if (z1 == 2)
  {
    #ifdef MYTRIM_ENABLED
        mooseError("alpha stopping not yet implemented!");
    #else
        std::cerr << "alpha stopping not yet implemented!\n";
        exit(1);
    #endif
  }
  else
  {
    yrmin = 0.13;
    vrmin = 1.0;
    v = std::sqrt(e / 25.0) / vfermi;

    if (v >= 1.0)
      vr = v * vfermi * (1.0 + 1.0 / (5.0 * v*v));
    else
      vr = (3.0 * vfermi / 4.0) * (1.0 + (2.0 * v*v / 3.0) - std::pow(v, 4.0) / 15.0);

    yr = std::max(yrmin, vr / std::pow(fz1, 0.6667));
    yr = std::max(yr, vrmin / std::pow(fz1, 0.6667));
    a = -0.803 * std::pow(yr, 0.3) + 1.3167 * std::pow(yr, 0.6) + 0.38157 * yr +  0.008983 * yr*yr;

    // ionization level of the ion at velocity yr
    q = std::min(1.0, std::max(0.0, 1.0 - std::exp(-std::min(a, 50.0))));

    b = (std::min(0.43, std::max(0.32, 0.12 + 0.025 * fz1))) / std::pow(fz1, 0.3333);
    l0 = (0.8 - q * std::min(1.2, 0.6 + fz1 / 30.0)) / std::pow(fz1, 0.3333);
    if (q < 0.2)
      l1 = 0.0;
    else if (q < std::max(0.0, 0.9 - 0.025 * fz1))
    {//210
      q1 = 0.2;
      l1 = b * (q - 0.2) / std::abs(std::max(0.0, 0.9 - 0.025 * fz1) - 0.2000001);
    }
    else if (q < std::max(0.0, 1.0 - 0.025 * std::min(16.0, fz1)))
      l1 = b;
    else
      l1 = b * (1.0 - q) / (0.025 * std::min(16.0, fz1));

    l = std::max(l1, l0 * lfctr);
    zeta = q + (1.0 / (2.0 * vfermi*vfermi)) * (1.0 - q) * std::log(1.0 + sqr(4.0 * l * vfermi / 1.919 ));

    // add z1^3 effect
    a = -sqr(7.6 - std::max(0.0, std::log(e)));
    zeta *= 1.0 + (1.0 / (fz1*fz1)) * (0.18 + 0.0015 * fz2) * std::exp(a);

    if (yr <= std::max(yrmin, vrmin / std::pow(fz1, 0.6667)))
    {
      // calculate velocity stopping for  yr < yrmin
      vrmin = std::max(vrmin, yrmin * std::pow(fz1, 0.6667));
      vmin = 0.5 * (vrmin + std::sqrt(std::max(0.0, vrmin*vrmin - 0.8 * vfermi*vfermi)));
      eee = 25.0 * vmin*vmin;
      sp = rpstop(z2, eee);

      if (z2 == 6 || ((z2 == 14 || z2 == 32) && z1 <= 19))
        power = 0.375;
      else
        power = 0.5;

      se = sp * sqr(zeta * fz1) * std::pow(e/eee, power);
    }
    else
    {
      sp = rpstop(z2, e);
      se = sp * sqr(zeta * fz1);
    }
  } // END: heavy-ions

  return se * 10.0;
}
