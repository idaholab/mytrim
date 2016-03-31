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
  const unsigned int end = _element.size();

  // get total stoichiometry
  for (unsigned int i = 0; i < end; ++i)
  {
    if (_element[i]->_t < 0.0)
      _element[i]->_t = 0.0;
    tt += _element[i]->_t;
  }

#ifdef MYTRIM_ENABLED
  if (tt == 0.0)
    mooseError("Stoichiometry invalid, all elements zero.");
#endif

  // normalize relative probabilities to 1
  for (unsigned int i = 0; i < end; ++i)
    _element[i]->_t /= tt;

  // average
  _am = 0.0;
  _az = 0.0;
  for (unsigned int i = 0; i < end; ++i)
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

  const unsigned int end = _element.size();
  for (unsigned int i = 0; i < end; ++i)
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
  const unsigned int end = _element.size();
  for (unsigned int i = 0; i < end; ++i)
    se += rstop(pka, _element[i]->_Z) * _element[i]->_t;

  return se * _arho;
}

Real
MaterialBase::rpstop(int z2p, Real e)
{
  Real pe, sl, sh, sp, velpwr;
  const int z2 = z2p - 1;

  // velocity proportional stopping below pe0
  const Real pe0 = 10.0; // [STO01210], 25.0 in original TRIM

  if (e > 1.0e4)
  {
    // high energy stopping
    const Real x = std::log(e) / e;
    const std::vector<Real> & ehigh = _simconf->scoef[z2-1].ehigh;

    sp = ehigh[0] + ehigh[1] * x + ehigh[2] * x*x + ehigh[3] / x;
  }
  else
  {
    pe = std::max(pe0, e);

    // pcoef indices are one less than in the fortran version!
    const std::vector<Real> & pcoef = _simconf->scoef[z2-1].pcoef;

    sl = pcoef[0] * std::pow(pe, pcoef[1]) + pcoef[2] * std::pow(pe, pcoef[3]);
    sh = pcoef[4] / std::pow(pe, pcoef[5]) * std::log(pcoef[6] / pe + pcoef[7] * pe);

    sp = sl * sh / (sl + sh);
    if (e <= pe0)
    {
      // velpwr is the power of velocity stopping below pe0
      if (z2p <= 6)
        velpwr = 0.35; // [STO01280], 0.25 in original TRIM
      else
        velpwr = 0.45;
      sp *= std::pow(e / pe0, velpwr);
    }
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
  unsigned int j;

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
    // Hydrogen electronic stopping powers [RST0640], pstop() in MCERD
    se = rpstop(z2, e);
  }
  else if (z1 == 2)
  {
    // Helium electronic stopping powers [RST0820], hestop() in MCERD
    const Real he0 = 1.0; // 10.0 in original TRIM

    Real he = std::max(he0, e);

    b = std::log(e);
    a = 0.2865 + 0.1266 * b - 0.001429 * b*b + 0.02402 * b*b*b - 0.1135 * std::pow(b, 4.0) + 0.001475 * std::pow(b, 5.0);

    Real heh = 1.0 - std::exp(-std::min(30.0, a));

    // add z1^3 effect to He/H stopping power ratio heh
    a = (1.0 + (0.007 + 0.00005 * z2) * std::exp(-std::pow((7.6 * std::max(0.0, std::log(he))), 2.0)));
    heh *= a * a;

    sp = rpstop(z2, he);
    se = sp * heh * 4.0;
    if (e <= he0)
      se *= std::sqrt(e / he0);
  }
  else
  {
    // Heavy ion electronic stopping powers [RST0990], histop() in MCERD
    yrmin = 0.13;
    vrmin = 1.0;
    const Real yrmin2 = vrmin / std::pow(fz1, 2.0/3.0);

    v = std::sqrt(e / 25.0) / vfermi;

    if (v >= 1.0)
      vr = v * vfermi * (1.0 + 1.0 / (5.0 * v*v));
    else
      vr = (3.0 * vfermi / 4.0) * (1.0 + (2.0 * v*v / 3.0) - std::pow(v, 4.0) / 15.0);

    yr = std::max(yrmin2, yrmin);
    yr = std::max(yr, vr / std::pow(fz1, 2.0/3.0));

    a = -0.803 * std::pow(yr, 0.3) + 1.3167 * std::pow(yr, 0.6) + 0.38157 * yr +  0.008983 * yr*yr;
    a = std::min(a, 50.0);

    // ionization level of the ion at velocity yr
    q = 1.0 - std::exp(a);
    q = q < 0.0 ? 0.0 : (q > 1.0 ? 1.0 : q);

    const std::vector<Real> & screen = _simconf->scoef[z1-1].screen;
    const std::vector<Real> & erange = _simconf->scoeflast.screen;

    for (j = 1; j <= 18 && q > erange[j]; ++j);
    j = j > 1 ? j - 1 : 1;
    j = j > 17 ? 17 : j;

    l0 = screen[j];
    // linear interpolation
    l1 = (q - erange[j]) * (screen[j+1] - screen[j]) /
                           (erange[j+1] - erange[j]);

    l = (l0 + l1) / std::pow(fz1, 1.0/3.0);

    zeta = q + (1.0 / (2.0 * vfermi*vfermi)) * (1.0 - q) * std::log(1.0 + sqr(4.0 * l * vfermi / 1.919));

    a = std::log(e);
    a = std::max(a, 0.0);

    zeta *= 1.0 + (1.0 / (fz1*fz1)) * (0.08 + 0.0015 * fz2) * std::exp(-sqr(7.6 - a));

    a = std::max(yrmin2, yrmin);

    if (yr <= a)
    {
      // calculate velocity stopping for  yr < yrmin
      vrmin = std::max(vrmin, yrmin * std::pow(fz1, 2.0/3.0));
      a = sqr(vmin) - 0.8 * sqr(vfermi);
      a = std::max(a, 0.0);

      vmin = 0.5 * (vrmin + std::sqrt(a));
      eee = 25.0 * vmin*vmin;
      sp = rpstop(z2, eee);
      Real eion = std::max(eee, 9999.0);

      const std::vector<Real> & vfcorr = _simconf->scoef[z2-1].fermicorr;
      const std::vector<Real> & vrange = _simconf->scoeflast.fermicorr;

      for (j = 1; j <= 13 && q > vrange[j]; ++j);
      j = j > 1 ? j - 1 : 1;
      j = j > 13 ? 13 : j;

      const Real vfcorr0 = vfcorr[j];
      // linear interpolation
      const Real vfcorr1 = (eion - vrange[j]) * (vfcorr[j+1] - vfcorr[j]) /
                                                (vrange[j+1] - vrange[j]);
      sp *= vfcorr0 + vfcorr1;

      power = 0.47;
      if (z1 == 3)
        power = 0.55;
      else if (z2 < 7)
        power = 0.375;
      else if (z1 < 18 && (z2 == 14 || z2 == 32))
        power = 0.375;

      se = sp * sqr(zeta * fz1) * std::pow(e/eee, power);
    }
    else
    {
      // sp = pstop(z2,E,scoef);
      // se = sp*intpow(zeta*z1,2);
      // eion = min(E,9999.0);
      // for(j=41;j<=53 && eion>=scoef[93][j];j++);
      // j--;
      // j = max(j,41);
      // j = min(j,53);
      //
      // vfcorr0 =scoef[z2][j];
      // vfcorr1 = (eion - scoef[93][j])*(scoef[z2][j+1] - scoef[z2][j])/
      //                                 (scoef[93][j+1] - scoef[93][j]);
      // se *= (vfcorr0 + vfcorr1);

      sp = rpstop(z2, e);
      se = sp * sqr(zeta * fz1);
    }
  } // END: heavy-ions

  return se * 10.0; // warum *10 ?!
}
