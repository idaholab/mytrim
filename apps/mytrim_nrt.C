#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <queue>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_wire.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"

#include "functions.h"

using namespace MyTRIM_NS;

#define _console std:: cout

SimconfType * _simconf;
MaterialBase * _material;

// some useful constants
static const Real math_pi = 3.14159265359;
static const Real electron_charge = 3.7941; // unit is sqrt(ev * A)

// composition of the problem
static const std::vector<std::string> names = {"C", "Si"};
static const std::vector<Real> A = {12.0, 28.0};
static const std::vector<Real> Z = {6.0, 14.0};
static const std::vector<Real> f = {0.5, 0.5};
static const std::vector<Real> displacement_threshold = {16.3, 92.6}; // from the paper

// simulation parameters
static const Real Emax = 1.0e3;
static const Real dE = 1.0;
static const bool use_power_law = true;
std::vector<Real> energy_nodes;
std::vector<std::vector<std::vector<Real>>> nu;

Real
stoppingPower(Real projectile_energy, unsigned int projectile_id)
{
  Real stopping_power = 0.0;
  Real weight = 0.0;

  // loop over all species and accumulate Se_composite using Bragg's Additivity rule
  IonBase * ion;
  for (unsigned int target_id = 0; target_id < A.size(); ++target_id)
  {
    // populate ion
    ion = new MyTRIM_NS::IonBase;
    ion->_gen = 0;
    ion->_tag = 0;
    ion->_Z = Z[projectile_id];
    ion->_m = A[projectile_id];
    ion->_E  = projectile_energy;

    ion->_dir(0) = 0.0;
    ion->_dir(1) = 1.0;
    ion->_dir(2) = 0.0;

    ion->_pos(0) = 0;
    ion->_pos(1) = 0;
    ion->_pos(2) = 0;

    int z2 = Z[target_id];
    Real s = _material->rstop(ion, z2);
    stopping_power += s * A[target_id] * f[target_id];
    weight += A[target_id] * f[target_id];
  }
  return stopping_power / weight;
}

Real
ft(Real t)
{
  Real lambda = 1.309;
  return lambda * pow(t, 1.0 / 6.0) * pow(1.0 + pow(2.0 * lambda * pow(t, 2.0 / 3.0), 2.0 / 3.0), -3.0 / 2.0);
}


Real
recoilCrossSection(Real projectile_energy, Real target_energy, unsigned int projectile_id, unsigned int target_id)
{
  Real Tm = 4.0 * A[target_id] * A[projectile_id] / sqr(A[target_id] + A[projectile_id]) * projectile_energy;
  if (target_energy > Tm)
    return 0;

  Real zt = Z[target_id];
  Real zp = Z[projectile_id];
  Real a = 0.4683 / sqrt(pow(zt, 2.0 / 3.0) + pow(zp, 2.0 / 3.0));
  Real kappa = projectile_energy * A[projectile_id] / A[target_id] * sqr(a / (2.0 * zt * zp * sqr(electron_charge)));
  Real t = target_energy * kappa;
  Real eta = projectile_energy * A[projectile_id] / A[target_id] * sqr(a / (2.0 * zt * zp * sqr(electron_charge)));
  return 0.5 * math_pi * sqr(a) * kappa * ft(t) / pow(t, 3.0 / 2.0);
}

void
getPieceWisePowerLaw(Real t, Real & alpha, Real & beta)
{
  if (t <= 1.0e-5)
  {
    alpha = 0;
    beta = 1;
  }
  else if (t > 1.0e-5 && t <= 5.0e-3)
  {
    alpha = -0.038867;
    beta = 0.6285;
  }
  else if (t > 5.0e-3 && t <= 1.0e-1)
  {
    alpha = -0.174285;
    beta = 0.30668;
  }
  else if (t > 1.0e-1 && t <= 1.0)
  {
    alpha = -0.354466;
    beta = 0.202539;
  }
  else if (t > 1.0    && t <= 10.0)
  {
    alpha = -0.50402;
    beta = 0.2025388;
  }
  else
  {
    alpha = -0.641031;
    beta = 0.3075704;
  }
}

Real
recoilCrossSectionPowerLaw(Real projectile_energy, Real target_energy, unsigned int projectile_id, unsigned int target_id)
{
  Real Tm = 4.0 * A[target_id] * A[projectile_id] / sqr(A[target_id] + A[projectile_id]) * projectile_energy;
  if (target_energy > Tm)
    return 0;

  Real zt = Z[target_id];
  Real zp = Z[projectile_id];
  Real At = A[target_id];
  Real Ap = A[projectile_id];
  Real lambda = 1.309;
  Real a0 = 0.4683;
  Real Eref = 2.0 * sqr(electron_charge) / a0;
  Real constant = 0.5 * math_pi * lambda * sqr(a0) / Eref;
  Real Lpt = Ap / At / (sqr(zt * zp) * (pow(zt, 2.0 / 3.0) + pow(zp, 2.0 / 3.0)));
  Real eps = projectile_energy / Eref;
  Real tau = target_energy / Eref;
  Real t =  Lpt * eps * tau;
  Real alpha, beta;
  getPieceWisePowerLaw(t, alpha, beta);
  return constant * Lpt / (pow(zt, 2.0 / 3.0) + pow(zp, 2.0 / 3.0)) * eps * beta * pow(t, alpha - 4.0 / 3.0);
}

void
checkData()
{
  std::vector<Real> energies = {0.15, 0.192, 0.24, 0.3, 0.384, 0.48, 0.6, 0.72, 0.84, 0.96, 1.08, 1.2,
                                1.5, 1.92, 2.4, 3.0, 3.84, 4.8, 6.0, 7.2, 8.4, 9.6, 10.8};
  for (auto & e : energies)
    _console << e << " MeV " << stoppingPower(e * 1.0e6, 0) / (10.0 * 19.94) << std::endl;

  _console << std::endl;
  for (auto & e : energies)
  {
    for (auto & t : energies)
      if (use_power_law)
        _console << recoilCrossSectionPowerLaw(e, t, 0, 0) << " ";
      else
        _console << recoilCrossSection(e, t, 0, 0) << " ";
    _console << std::endl;
  }
}

unsigned int
findIndex(Real energy)
{
  for (unsigned int ll = 0; ll < energy_nodes.size() - 1; ++ll)
    if (energy_nodes[ll + 1] > energy)
      return ll;
  return energy_nodes.size() - 1;
}

Real
getDeltaE(unsigned int l)
{
  if (l == 0)
    return 0.5 * energy_nodes[1]; // energy_nodes[0] = 0
  else if (l == energy_nodes.size() - 1)
    return 0.5 * (energy_nodes[l] - energy_nodes[l - 1]);
  return 0.5 * (energy_nodes[l + 1] - energy_nodes[l - 1]);
}

// computes the nuclear stopping power from sigma
Real
nuclearStopping(unsigned int projectile_id, unsigned int target_id, unsigned int l)
{
  Real integral = 0;
  Real Lij = 4.0 * A[target_id] * A[projectile_id] / sqr(A[target_id] + A[projectile_id]);

  // loop over all energies; E_0 = 0 < Ejd
  for (unsigned int ll = 1; ll < energy_nodes.size(); ++ll)
  {
    if (energy_nodes[ll] > Lij * energy_nodes[l])
      return integral;
    Real xs;
    if (use_power_law)
      xs = recoilCrossSectionPowerLaw(energy_nodes[l], energy_nodes[ll], projectile_id, target_id);
    else
      xs = recoilCrossSection(energy_nodes[l], energy_nodes[ll], projectile_id, target_id);
    integral += getDeltaE(l) * energy_nodes[ll] * xs;
  }
  return integral;
}

Real
integralI(unsigned int projectile_id, unsigned int target_id, unsigned int l)
{
  Real integral = 0;
  Real Lij = 4.0 * A[target_id] * A[projectile_id] / sqr(A[target_id] + A[projectile_id]);

  // check easy way out
  if (energy_nodes[l] < displacement_threshold[target_id])
    return 0;

  // loop over all energies; E_0 = 0 < Ejd
  for (unsigned int ll = 1; ll < energy_nodes.size(); ++ll)
  {
    if (energy_nodes[ll] > Lij * energy_nodes[l])
      return integral;

    Real xs;
    if (use_power_law)
      xs = recoilCrossSectionPowerLaw(energy_nodes[l], energy_nodes[ll], projectile_id, target_id);
    else
      xs = recoilCrossSection(energy_nodes[l], energy_nodes[ll], projectile_id, target_id);

    if (energy_nodes[ll] > displacement_threshold[target_id])
      integral += getDeltaE(ll) * xs;
  }
  return integral;
}

Real
integralII(unsigned int projectile_id, unsigned int target_id, unsigned int k, unsigned int l)
{
  Real integral = 0;
  Real Lik = 4.0 * A[k] * A[projectile_id] / sqr(A[k] + A[projectile_id]);

  // check easy way out
  if (energy_nodes[l] < displacement_threshold[target_id])
    return 0;

  // loop over all energies; E_0 = 0 implies nu_ij = 0
  for (unsigned int ll = 1; ll < energy_nodes.size(); ++ll)
  {
    if (energy_nodes[ll] > Lik * energy_nodes[l])
      return integral;

    Real xs;
    if (use_power_law)
      xs = recoilCrossSectionPowerLaw(energy_nodes[l], energy_nodes[ll], projectile_id, k);
    else
      xs = recoilCrossSection(energy_nodes[l], energy_nodes[ll], projectile_id, k);

    if (ll == l)
      // in this case we don't know the nu at energy l since we are just computing it.
      // constant approximation consistent with forward Euler tells us to use value at l-1
      integral += getDeltaE(ll) * nu[k][target_id][ll - 1] * xs;
    else if (ll != l && energy_nodes[ll] > displacement_threshold[target_id])
      integral += getDeltaE(ll) * nu[k][target_id][ll] * xs;
  }
  return integral;
}

Real
integralIII(unsigned int projectile_id, unsigned int target_id, unsigned int k, unsigned int l)
{
  Real integral = 0;
  Real Lik = 4.0 * A[k] * A[projectile_id] / sqr(A[k] + A[projectile_id]);
  Real Lij = 4.0 * A[target_id] * A[projectile_id] / sqr(A[target_id] + A[projectile_id]);

  // we cannot allow T = 0 because the scattering XS blows up
  // in this case no energy is lost so this is an uninteresting event (is that true)?
  for (unsigned int ll = 1; ll < energy_nodes.size(); ++ll)
  {
    if (energy_nodes[ll] > Lik * energy_nodes[l])
      return integral;

    Real Ediff = energy_nodes[l] - energy_nodes[ll];

    // this takes care of the heavyside step function!
    if (Ediff - displacement_threshold[target_id] / Lij < 0.0)
      continue;

    Real xs;
    if (use_power_law)
      xs = recoilCrossSectionPowerLaw(energy_nodes[l], energy_nodes[ll], projectile_id, k);
    else
      xs = recoilCrossSection(energy_nodes[l], energy_nodes[ll], projectile_id, k);

    // find interpolated nu => simply take the value from the left
    unsigned int index = findIndex(Ediff);
    Real interpolated_nu = nu[projectile_id][target_id][index];

    integral += getDeltaE(ll) * interpolated_nu * xs;
  }
  return integral;
}

Real
integralIV(unsigned int projectile_id, unsigned int target_id, unsigned int k, unsigned int l)
{
  Real integral = 0;
  Real Lik = 4.0 * A[k] * A[projectile_id] / sqr(A[k] + A[projectile_id]);

  // loop over all energies; E_0 = 0 < Ejd
  for (unsigned int ll = 1; ll < energy_nodes.size(); ++ll)
  {
    if (energy_nodes[ll] > Lik * energy_nodes[l])
      return integral;

      Real xs;
      if (use_power_law)
        xs = recoilCrossSectionPowerLaw(energy_nodes[l], energy_nodes[ll], projectile_id, k);
      else
        xs = recoilCrossSection(energy_nodes[l], energy_nodes[ll], projectile_id, k);

    integral += getDeltaE(ll) * nu[projectile_id][target_id][l - 1] * xs;
  }
  return integral;
}

Real
nrtRightHandSide(unsigned int projectile_id, unsigned int target_id, unsigned int l)
{
  Real fijl = 0;
  // loop over k as in Eq. (7)
  for (unsigned int k = 0; k < A.size(); ++k)
  {
    // NOTE: j == k!
    // integralI = integral_Ejd^{Lambda_ij * E_l} dT sigma_ij(E,T)
    if (target_id == k)
      fijl += f[k] * integralI(projectile_id, target_id, l);

    // integralII_ijkl = integral_{Ekd}^{Lambda_ik E_l} dT nu_kj * sigma_ik(E_l, T)
    fijl += f[k] * integralII(projectile_id, target_id, k, l);

    // integralIII_ijkl = integral_{0}^{Lambda_ik} Gamma(E-T-Ejd/Lambda_ij) nu_ij(E_l - T) sigma_ik(E_l, T) dT
    fijl += f[k] * integralIII(projectile_id, target_id, k, l);

    // integralIV_ijkl = integral_{0}^{Lambda_ij * E_l} sigma_ik(E_l, T) * vij(E) dT
    // Note: it's subtracted!
    fijl -= f[k] * integralIV(projectile_id, target_id, k, l);
  }
  return fijl;
}

void
takeNRTStep(unsigned int l)
{
  for (unsigned int projectile_id = 0; projectile_id < A.size(); ++projectile_id)
    for (unsigned int target_id = 0; target_id < A.size(); ++target_id)
    {
      Real stopping_power = stoppingPower(energy_nodes[l], projectile_id);
      Real deltaE = energy_nodes[l] - energy_nodes[l - 1];
      nu[projectile_id][target_id][l] = nu[projectile_id][target_id][l - 1] + deltaE / stopping_power
                                      * nrtRightHandSide(projectile_id, target_id, l);
    }
}

int
main()
{
  unsigned int seed = 1;

  // initialize global parameter structure and read data tables from file
  _simconf = new SimconfType(seed);

  _material = new MaterialBase(_simconf, 1.0);

  // check out data
  //checkData();

  // set up energy_nodes: assume equal spacing
  Real Ec = 0.0;
  while (Ec < Emax)
  {
    energy_nodes.push_back(Ec);
    Ec += dE;
  }

  // set up nu
  nu.resize(A.size());
  for (unsigned int j = 0; j < A.size(); ++j)
  {
    nu[j].resize(A.size());
    for (unsigned int k = 0; k < A.size(); ++k)
      nu[j][k].resize(energy_nodes.size());
  }

  // step through all times and compute nu_ij for that time
  // Note: no need to work on l = 0 since nu_ij = 0
  for (unsigned int l = 1; l < energy_nodes.size(); ++l)
    takeNRTStep(l);

  for (unsigned int l = 0; l < energy_nodes.size(); ++l)
  {
    for (unsigned int i = 0; i < A.size(); ++i)
      for (unsigned int j = 0; j < A.size(); ++j)
        _console << energy_nodes[l] << " " << nu[i][j][l] << " ";
    _console << std::endl;
  }

  delete _simconf;
  delete _material;
}
