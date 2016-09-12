#include "material.h"
#include "simconf.h"

#include "functions.h"

#include <cmath>
#include <iostream>

using namespace MyTRIM_NS;

MaterialBase::MaterialBase(SimconfType * simconf, Real rho) :
    _rho(rho),
    _tag(-1),
    _dirty(true),
    _simconf(simconf)
{
}

void
MaterialBase::prepare()
{
  Real tt = 0.0;
  const unsigned int end = _element.size();

  // get total stoichiometry
  for (unsigned int i = 0; i < end; ++i)
  {
    if (_element[i]._t < 0.0)
      _element[i]._t = 0.0;
    tt += _element[i]._t;
  }

#ifdef MYTRIM_ENABLED
  if (tt == 0.0)
    mooseError("Stoichiometry invalid, all elements zero.");
#endif

  // normalize relative probabilities to 1
  for (unsigned int i = 0; i < end; ++i)
    _element[i]._t /= tt;

  // average
  _am = 0.0;
  _az = 0.0;
  for (unsigned int i = 0; i < end; ++i)
  {
    _am += _element[i]._m * _element[i]._t;
    _az += Real(_element[i]._Z) * _element[i]._t;
  }

  _arho = _rho * 0.6022 / _am; //[TRI00310] atoms/Ang^3
}

// make sure layers are prepare'd first!
void
MaterialBase::average(const IonBase * pka)
{
  mu = pka->_m / _am;
  const Real fZ = Real(pka->_Z);

  // universal or firsov screening length
  const Real fZ023 = std::pow(fZ, 0.23);
  a = .5292 * .8853 / (fZ023 + std::pow(_az, 0.23));
  //a = .5292 * .8853 / std::pow(pow(Real(pka._Z), 0.5) + std::pow(_az, 0.5), 2.0/3.0);

  // mean flight path0
  f = a * _am / (_az * fZ * 14.4 * (pka->_m + _am));
  //eps0 = e0 * f;
  epsdg = _simconf->tmin * f * Utility::pow<2>(1.0 + mu) / (4.0 * mu);

  // fd and kd determine how much recoil energy goes into el. loss and vaccancies
  fd = std::pow(0.01 * _az, -7.0 / 3.0);
  kd = std::pow(0.1334 * _az, 2.0 / 3.0) / std::sqrt(_am);

  const unsigned int end = _element.size();
  for (unsigned int i = 0; i < end; ++i)
  {
    _element[i].my = pka->_m / _element[i]._m;
    _element[i].ec = 4.0 * _element[i].my / Utility::pow<2>(1.0 + _element[i].my);
    _element[i].ai = .5292 * .8853 / (fZ023 + std::pow(_element[i]._Z, 0.23));
    //ai = .5292 * .8853 / std::pow(pow(Real(pka._Z), 0.5) + std::pow(_element[i].z, 0.5), 2.0/3.0);
    _element[i].fi = _element[i].ai * _element[i]._m /
                     (fZ * Real(_element[i]._Z) * 14.4 * (pka->_m + _element[i]._m));
  }

  _dirty = false;
}

// make sure layers are prepare'd and averaged first!
Real
MaterialBase::getrstop(const IonBase * pka)
{
  Real se = 0.0;
  const unsigned int end = _element.size();
  for (unsigned int i = 0; i < end; ++i)
    se += rstop(pka, _element[i]._Z) * _element[i]._t;

  return se * _arho;
}

Real
MaterialBase::rpstop(int z2p, Real e)
{
  Real pe, sl, sh, sp, velpwr;
  const int z2 = z2p - 1;
  // velocity proportional stopping below pe0
  const Real pe0 = 25.0;
  pe = std::max(pe0, e);

  // pcoef indices are one less than in the fortran version!
  sl = (_simconf->scoef[z2].pcoef[0] * std::pow(pe, _simconf->scoef[z2].pcoef[1])) +
       (_simconf->scoef[z2].pcoef[2] * std::pow(pe, _simconf->scoef[z2].pcoef[3]));
  sh = _simconf->scoef[z2].pcoef[4] / std::pow(pe, _simconf->scoef[z2].pcoef[5]) *
       std::log(_simconf->scoef[z2].pcoef[6] / pe + _simconf->scoef[z2].pcoef[7] * pe);
  sp = sl * sh / (sl + sh);
  if (e <= pe0)
  {
    // velpwr is the power of velocity stopping below pe0
    if (z2p <= 6)
      velpwr = 0.25;
    else
      velpwr = 0.45;
    sp *= std::pow(e / pe0, velpwr);
  }
  return sp;
}

Real
MaterialBase::rstop(const IonBase * ion, int z2)
{
  Real e, vrmin, yrmin, v, vr, yr, vmin, m1;
  Real a, b, q, /*q1,*/ l, l0, l1;
  Real zeta;
  const int z1 = ion->_Z;
  const Real fz1 = Real(z1);
  const Real fz2 = Real(z2);
  Real eee, sp, power;
  Real se;

  // scoeff
  const Real lfctr = _simconf->scoef[z1-1].lfctr;
  const Real mm1 = _simconf->scoef[z1-1].mm1;
  const Real vfermi = _simconf->scoef[z2-1].vfermi;
  //Real atrho = _simconf->scoef[z2-1].atrho;

  if (ion->_m == 0.0)
    m1 = mm1;
  else
    m1 = ion->_m;

  // we store ion energy in eV but ee is needed in keV
  const Real ee = 0.001 * ion->_E;
  e = ee / m1;

  if (z1 == 1)
  {
    // Hydrogen electronic stopping powers [RST0640]
    se = rpstop(z2, e);
  }
  else if (z1 == 2)
  {
    // Helium electronic stopping powers [RST0820]
    const Real he0 = 1.0;
    Real he = std::max(he0, e);

    b = std::log(he);
    const Real b2 = b*b;
    const Real b4 = b2*b2;
    a = 0.2865 + 0.1266 * b - 0.001429 * b2 + 0.02402 * b*b2 - 0.01135 * b4 + 0.001475 * b4*b;
    Real heh = 1.0 - std::exp(-std::min(30.0, a));

    he = std::max(he, 1.0);
    a = 1.0 + (0.007 + 0.00005 * fz2) * std::exp(-Utility::pow<2>(7.6 - std::log(he)));
    heh *= a * a;

    sp = rpstop(z2, he);
    se = sp * heh * 4.0;
    if (e <= he0)
      se *= std::sqrt(e / he0);
  }
  else
  {
    // Heavy ion electronic stopping powers [RST0990]
    yrmin = 0.13;
    vrmin = 1.0;

    v = std::sqrt(e / 25.0) / vfermi;
    const Real v2 = v*v;

    if (v >= 1.0)
      vr = v * vfermi * (1.0 + 1.0 / (5.0 * v2));
    else
      vr = (3.0 * vfermi / 4.0) * (1.0 + (2.0 * v2 / 3.0) - v2*v2 / 15.0);

    const Real cbrt_fz1 = std::cbrt(fz1);
    const Real cbrt2_fz1 = cbrt_fz1 * cbrt_fz1;
    yr = std::max(yrmin, vr / cbrt2_fz1);
    yr = std::max(yr, vrmin / cbrt2_fz1);
    const Real yr03 = std::pow(yr, 0.3);
    a = -0.803 * yr03 + 1.3167 * yr03*yr03 + 0.38157 * yr +  0.008983 * yr*yr;

    // ionization level of the ion at velocity yr
    q = std::min(1.0, std::max(0.0, 1.0 - std::exp(-std::min(a, 50.0))));

    b = (std::min(0.43, std::max(0.32, 0.12 + 0.025 * fz1))) / cbrt_fz1;
    l0 = (0.8 - q * std::min(1.2, 0.6 + fz1 / 30.0)) / cbrt_fz1;
    if (q < 0.2)
      l1 = 0.0;
    else if (q < std::max(0.0, 0.9 - 0.025 * fz1))
    {//210
      // q1 = 0.2; in the original code, but never used
      l1 = b * (q - 0.2) / std::abs(std::max(0.0, 0.9 - 0.025 * fz1) - 0.2000001);
    }
    else if (q < std::max(0.0, 1.0 - 0.025 * std::min(16.0, fz1)))
      l1 = b;
    else
      l1 = b * (1.0 - q) / (0.025 * std::min(16.0, fz1));

    l = std::max(l1, l0 * lfctr);
    zeta = q + (1.0 / (2.0 * vfermi*vfermi)) * (1.0 - q) * std::log(1.0 + Utility::pow<2>(4.0 * l * vfermi / 1.919 ));

    // add z1^3 effect
    a = -Utility::pow<2>(7.6 - std::max(0.0, std::log(e)));
    zeta *= 1.0 + (1.0 / (fz1*fz1)) * (0.18 + 0.0015 * fz2) * std::exp(a);

    if (yr <= std::max(yrmin, vrmin / cbrt2_fz1))
    {
      // calculate velocity stopping for  yr < yrmin
      vrmin = std::max(vrmin, yrmin * cbrt2_fz1);
      vmin = 0.5 * (vrmin + std::sqrt(std::max(0.0, vrmin*vrmin - 0.8 * vfermi*vfermi)));
      eee = 25.0 * vmin*vmin;
      sp = rpstop(z2, eee);

      if (z2 == 6 || ((z2 == 14 || z2 == 32) && z1 <= 19))
        power = 0.375;
      else
        power = 0.5;

      se = sp * Utility::pow<2>(zeta * fz1) * std::pow(e/eee, power);
    }
    else
    {
      sp = rpstop(z2, e);
      se = sp * Utility::pow<2>(zeta * fz1);
    }
  } // END: heavy-ions

  return se * 10.0;
}
