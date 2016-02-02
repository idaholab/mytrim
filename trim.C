#include "trim.h"

#include "functions.h"
#include <iostream>

using namespace MyTRIM_NS;

//#define RANGECORRECT2

//
// all energies are in eV
//

// does a single ion cascade
void trimBase::trim(ionBase *pka_, std::queue<ionBase*> &recoils)
{
  // simconf should already be initialized
  pka = pka_;

  // make recoil queue available in overloadable functions
  recoil_queue_ptr = &recoils;

  //e = pka.e;
  Real pl = 0.0;
  Real max = 0.0;
  //Real e0kev = pka->e / 1000.0;
  int ic = 0;
  unsigned int nn; //, ie;
  Real r1, r2, hh;
  Real eps, eeg, ls, p, b, r, see, dee;
  Real s2, c2, ct, st;
  Real rr, ex1, ex2, ex3, ex4, v ,v1;
  Real fr, fr1, q, roc, sqe;
  Real cc, aa, ff, co, delta;
  Real den;
  Real rdir[3], perp[3], norm, psi;

  Real p1, p2;
  //Real range;

  // generate random number for use in the first loop iteration only!
  r1 = dr250();

  do // cycle for each collision
  {
    // increase loop counter
    ic++;

    // which material is the ion currently in?
    material = sample->lookupMaterial(pka->pos);
    if (material==0) break; // TODO: add flight through vacuum

    // normalize direction vector
    v_norm(pka->dir);

    // setup max. impact parameter
    eps =  pka->e * material->f;
    eeg = std::sqrt(eps * material->epsdg); // [TRI02450]
    material->pmax = material->a / (eeg + std::sqrt(eeg) + 0.125 * std::pow(eeg, 0.1));

    ls = 1.0 / (M_PI * std::pow(material->pmax, 2.0) * material->arho);
    if (ic==1) ls = r1 * fmin(ls, simconf->cw);

    // correct for maximum available range in current material by increasing maximum impact parameter
    #ifdef RANGECORRECT
      range = sample->rangeMaterial(pka->pos, pka->dir);
      if (range<ls)
      {
        /* std::cout << "range=" << range << " ls=" << ls
              << " pos[0]=" << pka->pos[0] << " dir[0]=" << pka->dir[0] << std::endl;
        std::cout << "CC " << pka->pos[0] << ' ' << pka->pos[1] << std::endl;
        std::cout << "CC " << pka->pos[0] + pka->dir[0] * range << ' ' << pka->pos[1] + pka->dir[1] * range << std::endl;
        std::cout << "CC " << std::endl;*/

        ls = range;

        // correct pmax to correspond with new ls
        material->pmax = 1.0 / std::sqrt( M_PI * ls * material->arho);
      }
    #endif

    // correct for maximum available range in current material by dropping recoils randomly (faster)
    #ifdef RANGECORRECT2
      range = sample->rangeMaterial(pka->pos, pka->dir);
      if (range < ls)
      {
        // skip this recoil, just advance the ion
        if (range/ls < dr250())
        {
          // electronic stopping
          pka->e -= range * material->getrstop(pka);

          // free flight
          for (int i = 0; i < 3; i++)
            pka->pos[i] += pka->dir[i] * range;

          // start over
          continue;
        }
        ls = range;
      }
    #endif

    // advance clock pathlength/velocity
    pka->t += 10.1811859 * (ls - simconf->tau) / std::sqrt(2.0 * pka->e / pka->_m);

    // time in fs! m in u, l in Ang, e in eV
    // 1000g/kg, 6.022e23/mol, 1.602e-19J/eV, 1e5m/s=1Ang/fs 1.0/0.09822038
    //printf("se %d  %f [eV]  %f [keV/nm]  %f [nm]\n", pka->id, pka->e, see/100.0, pl/10.0);

    // choose impact parameter
    r2 = dr250();
    p = material->pmax * std::sqrt(r2);

    // which atom in the material will be hit
    hh = dr250(); // selects element inside material to scatter from
    for (nn = 0; nn < material->element.size(); ++nn)
    {
      hh -= material->element[nn]->t;
      if (hh <= 0) break;
    }
    element = material->getElement(nn);

    // epsilon and reduced impact parameter b
    eps = element->fi * pka->e;
    b = p / element->ai;

    ////ie = int(pka.e / e0kev - 0.5); // was +0.5 for fortran indices
    //ie = int(pka.e / material->semax - 0.5); // was +0.5 for fortran indices
    //see = material->se[ie];

    see = material->getrstop(pka);
    //if (pka.e < e0kev) see = material->se[0] * std::sqrt(pka.e / e0kev);
    dee = ls * see;

    if (eps>10.0)
    {
      // use rutherford scattering
      s2 = 1.0 / (1.0 + (1.0 + b * (1.0 + b)) * std::pow(2.0 * eps * b , 2.0));
      c2 = 1.0 - s2;
      ct = 2.0 * c2 - 1.0;
      st = std::sqrt(1.0 - ct*ct);
    }
    else
    {
      // first guess at ion c.p.a. [TRI02780]
      r = b;
      rr = -2.7 * logf(eps * b);
      if (rr >= b)
      {
        r = rr;
        rr = -2.7 * logf(eps * rr);
        if (rr >= b) r = rr;
      }

      do
      {
        // universal potential
        ex1 = 0.18175 * exp(-3.1998 * r);
        ex2 = 0.50986 * exp(-0.94229 * r);
        ex3 = 0.28022 * exp(-0.4029 * r);
        ex4 = 0.028171 * exp(-0.20162 * r);
        v = (ex1 + ex2 + ex3 + ex4) / r;
        v1 = -(v + 3.1998 *ex1 + 0.94229 * ex2 + 0.4029 * ex3 + 0.20162 * ex4) / r;

        fr = b*b / r + v * r / eps -r;
        fr1 = - b*b / (r*r) + (v + v1 * r) / eps - 1.0;
        q = fr / fr1;
        r -= q;
      } while (std::abs(q/r) > 0.001); // [TRI03110]

      roc = -2.0 * (eps - v) / v1;
      sqe = std::sqrt(eps);

      // 5-parameter magic scattering calculation (universal pot.)
      cc = (0.011615 + sqe) / (0.0071222 + sqe);               // 2-87 beta
      aa = 2.0 * eps * (1.0 + (0.99229 / sqe)) * std::pow(b, cc); // 2-87 A
      ff = (std::sqrt(aa*aa + 1.0) - aa) * ((9.3066 + eps) / (14.813 + eps));

      delta = (r - b) * aa * ff / (ff + 1.0);
      co = (b + delta + roc) / (r + roc);
      c2 = co*co;
      s2 = 1.0 - c2;
      //printf("nonrf\n");
      ct = 2.0 * c2 - 1.0;
      st = std::sqrt(1.0 - ct*ct);
    } // end non-rutherford scattering

    // energy transferred to recoil atom
    den = element->ec * s2 * pka->e;

    if (dee > pka->e) {
      // avoid getting negative energies
      dee = pka->e;

      // sanity check
      if (den > 100.0)
        std::cerr << " electronic energy loss stopped the ion. Broken recoil!!\n";
    }

    // electronic energy loss
    pka->e -= dee;
    simconf->EelTotal += dee;

    // momentum transfer
    p1 = std::sqrt(2.0 * pka->_m * pka->e); // momentum before collision
    if (den > pka->e) den = pka->e; // avoid negative energy
    pka->e -= den;
    p2 = std::sqrt(2.0 * pka->_m * pka->e); // momentum after collision

    // track maximum electronic energy loss TODO: might want to track max(see)!
    if (dee>max) max = dee;

    // total path lenght
    pl += ls - simconf->tau;

    // find new position, save old direction to recoil
    recoil = pka->spawnRecoil();
    for (int i = 0; i < 3; i++)
    {
      // used to assign the new position to the recoil, but
      // we have to make sure the recoil starts in the appropriate material!
      pka->pos[i] += pka->dir[i] * (ls - simconf->tau);
      recoil->dir[i] = pka->dir[i] * p1;
    }
    recoil->e = den;

    // recoil loses the lattice binding energy
    recoil->e -= element->Elbind;
    recoil->_m = element->m;
    recoil->_Z = element->z;

    // create a random vector perpendicular to pka.dir
    // there is a cleverer way by using the azimuthal angle of scatter...
    do
    {
      for (int i = 0; i < 3; i++) rdir[i] = dr250() - 0.5;
      v_cross(pka->dir, rdir, perp);
      norm = std::sqrt(v_dot(perp, perp));
    }
    while (norm == 0.0);
    v_scale(perp, 1.0 / norm);

    psi = atan(st / (ct + element->my));
    v_scale(pka->dir, std::cos(psi));

    // calculate new direction, subtract from old dir (stored in recoil)
    for (int i = 0; i < 3; i++)
    {
      pka->dir[i] += perp[i] * std::sin(psi);
      recoil->dir[i] -= pka->dir[i] * p2;
    }

    // end cascade if a CUT boundary is crossed
    for (int i = 0; i < 3; i++) {
      if (sample->bc[i]==sampleBase::CUT &&
           (pka->pos[i]>sample->w[i] || pka->pos[i]<0.0)) {
        pka->state = ionBase::LOST;
        break;
      }
    }

    //
    // decide on the fate of recoil and pka
    //
    if (pka->state != ionBase::LOST) {
      if (recoil->e > element->Edisp) {
        // non-physics based descision on recoil following
        if (followRecoil()) {
          v_norm(recoil->dir);
          recoil->tag = material->tag;
          recoil->id  = simconf->id++;

          // queue recoil for processing
          recoils.push(recoil);
          if (simconf->fullTraj)
            std::cout << "spawn " << recoil->id << ' ' << pka->id << std::endl;
        } else {
          // this recoil could have left its lattice site, but we chose
          // not to follow it (simulation of PKAs only)
          recoil->id = ionBase::DELETE;
        }

        // will the knock-on get trapped at the recoil atom site?
        // (TODO: make sure that pka->ef < element->Edisp for all elements!)
        // did we create a vacancy by knocking out the recoil atom?
        if (pka->e > element->Edisp) {
          // yes, because the knock-on can escape, too!
          vacancyCreation();
        } else {
          // nope, the pka gets stuck at that site as...
          if (pka->_Z == element->z)
            pka->state = ionBase::REPLACEMENT;
          else
            pka->state = ionBase::SUBSTITUTIONAL;
        }
      } else {
        // this recoil will not leave its lattice site
        dissipateRecoilEnergy();
        recoil->id  = ionBase::DELETE;;

        // if the PKA has no energy left, put it to rest here as an interstitial
        if (pka->e < pka->ef) {
          pka->state = ionBase::INTERSTITIAL;
        }
      }
    }

    // delete recoil if it was not queued
    if (recoil->id  == ionBase::DELETE) delete recoil;

    // act on the pka state change
    checkPKAState();

    // output the full trajectory (state is not output by the ion object)
    if (simconf->fullTraj)
      std::cout << pka->state << ' ' << *pka << std::endl;

  } while (pka->state == ionBase::MOVING);
}

void trimBase::vacancyCreation()
{
  simconf->vacancies_created++;

  /*
  // modified Kinchin-Pease
  if (recoil->gen == 1)
  {
    // calculate modified kinchin pease data http://www.iue.tuwien.ac.at/phd/hoessinger/node47.html
    Real ed = 0.0115 * std::pow(material->az, -7.0/3.0) * recoil->e;
    Real g = 3.4008 * std::pow(ed, 1.0/6.0) + 0.40244 * std::pow(ed, 3.0/4.0) + ed;
    Real kd = 0.1337 * std::pow(material->az, 2.0/3.0) / std::pow(material->am, 0.5); //Z,M
    Real Ev = recoil->e / (1.0 + kd * g);
    simconf->KP_vacancies += 0.8 * Ev / (2.0 * element->Edisp);
    // should be something like material->Edisp (average?)
  }
  */
}


/*
materialBase* sampleType::lookupLayer(const Real* pos)
{
  Real dif[3];

  dif[0] = pos[0] - 100.0;
  dif[1] = pos[1];
  dif[2] = pos[2];
  Real r2 = v_dot(dif, dif);
  if (r2 < 2500.0) // r<50.0
    return material[1];
  else
    return material[0];
}
*/
