#include <math.h>
#include <stdio.h>

#include "trim.h"
#include "simconf.h"
#include <r250.h>

#include "functions.h"

//
// all energies are in eV
//

// does a single ion cascade
void trimBase::trim( ionBase *pka_, queue<ionBase*> &recoils )
{
  // simconf should already be initialized
  pka = pka_;

  //e = pka.e;
  float pl = 0.0;
  float max = 0.0;
  float e0kev = pka->e / 1000.0;
  int ic = 0;
  int nn, ie;
  double r1, r2, hh;
  float eps, eeg, ls, p, b, r, see, dee;
  float s2, c2, ct, st;
  float rr, ex1, ex2, ex3, ex4, v ,v1;
  float fr, fr1, q, roc, sqe;
  float cc, aa, ff, co, delta;
  float den;
  float rdir[3], perp[3], norm, psi;

  float p1, p2;
  float range;
  bool terminate;
//if( simconf->fullTraj )
  r1 = dr250();

  do // cycle for each collision
  {
    r2 = dr250();
    hh = dr250(); // selects element inside material to scatter from

    ic++;
    material = sample->lookupMaterial( pka->pos );
    if( material == 0 ) break;

    eps =  pka->e * material->f;
    // setup max. impact parameter
    eeg = sqrtf( eps * material->epsdg );
    material->pmax = material->a / ( eeg + sqrtf( eeg ) + 0.125 * pow( eeg, 0.1 ) );
    ls = 1.0 / ( M_PI * pow( material->pmax, 2.0 ) * material->arho );

    if( ic == 1 ) ls = r1 * fmin( ls, simconf->cw );

    // correct for maximum available range in current material
    range = sample->rangeMaterial( pka->pos, pka->dir );
    if( range < ls )
    {
      ls = range;
      // correct pmax to correspond with new ls
      material->pmax = 1.0 / sqrtf(  M_PI * ls * material->arho );
    }

    p = material->pmax * sqrtf( r2 );

    // if(hh>1.0 || r2>1.0 ) { fprintf(stderr,"weird random number!\n");}

    // which atom in the material will be hit
    for( nn = 0; nn < material->element.size(); nn++ )
    {
      hh -= material->element[nn]->t;
      if( hh <= 0 ) break;
    }
    // epsilon and reduced impact parameter b
    eps = material->element[nn]->fi * pka->e;
    b = p / material->element[nn]->ai;

    ////ie = int( pka.e / e0kev - 0.5 ); // was +0.5 for fortran indices
    //ie = int( pka.e / material->semax - 0.5 ); // was +0.5 for fortran indices
    //see = material->se[ie];

    see = material->getrstop( pka );
    //if( pka.e < e0kev ) see = material->se[0] * sqrtf( pka.e / e0kev );
    dee = ls * see;

    if( eps > 10.0 )
    {
      // use rutherford scattering
      s2 = 1.0 / ( 1.0 + ( 1.0 + b * ( 1.0 + b ) ) * pow( 2.0 * eps * b , 2.0 ) );
      c2 = 1.0 - s2;
      ct = 2.0 * c2 - 1.0;
      st = sqrtf( 1.0 - ct*ct );
    }
    else
    {
      // first gues at ion c.p.a. [TRI02780]
      r = b;
      rr = -2.7 * logf( eps * b );
      if( rr >= b )
      {
        r = rr;
        rr = -2.7 * logf( eps * rr );
        if( rr >= b ) r = rr;
      }

      do
      {
        // universal potential
        ex1 = 0.18175 * exp( -3.1998 * r );
        ex2 = 0.50986 * exp( -0.94229 * r );
        ex3 = 0.28022 * exp( -0.4029 * r );
        ex4 = 0.028171 * exp( -0.20162 * r );
        v = ( ex1 + ex2 + ex3 + ex4 ) / r;
        v1 = -( v + 3.1998 *ex1 + 0.94229 * ex2 + 0.4029 * ex3 + 0.20162 * ex4 ) / r;

        fr = b*b / r + v * r / eps -r;
        fr1 = - b*b / ( r*r ) + ( v + v1 * r ) / eps - 1.0;
        q = fr / fr1;
        r -= q;
      } while( fabs( q / r ) > 0.001 ); // [TRI03110]

      roc = -2.0 * ( eps - v ) / v1;
      sqe = sqrtf( eps );

      // 5-parameter magic scattering calculation (universal pot.)
      cc = ( 0.011615 + sqe ) / ( 0.0071222 + sqe );
      aa = 2.0 * eps * ( 1.0 + ( 0.99229 / sqe ) ) * pow( b, cc ); 
      ff = ( sqrtf( aa*aa + 1.0 ) - aa ) * ( ( 9.3066 + eps ) / ( 14.813 + eps ) );

      delta = ( r - b ) * aa * ff / ( ff + 1.0 );
      co = ( b + delta + roc ) / ( r + roc );
      c2 = co*co;
      s2 = 1.0 - c2;
//printf("nonrf\n");
      ct = 2.0 * c2 - 1.0;
      st = sqrtf( 1.0 - ct*ct );
    } // end non-rutherford scattering

    // energy transferred to recoil atom
    den = material->element[nn]->ec * s2 * pka->e;

    // advance clock pathlength/velocity
    pka->t += 10.1811859 * ( ls - simconf->tau ) / sqrt( 2.0 * pka->e / pka->m1 );
    // time in fs! m in u, l in Ang, e in eV 
    // 1000g/kg, 6.022e23/mol, 1.602e-19J/eV, 1e5m/s=1Ang/fs 1.0/0.09822038
    //printf( "se %d  %f [eV]  %f [keV/nm]  %f [nm]\n", pka->id, pka->e, see/100.0, pl/10.0 );

    pka->e -= dee; // electronic energy loss
    if( pka->e < 0.0 && den > 100.0 ) 
      fprintf( stderr, " electronic energy loss stopped the ion. Broken recoil!!\n" );

    p1 = sqrtf( 2.0 * pka->m1 * pka->e ); // momentum before collision
    pka->e -= den; if( pka->e < 0.0 ) pka->e = 0.0;
    p2 = sqrtf( 2.0 * pka->m1 * pka->e ); // momentum after collision 

    if( dee > max ) max = dee;

    // total path lenght
    pl += ls - simconf->tau;

    // find new position, save old direction to recoil
    v_norm( pka->dir );
    recoil = new ionBase;
    for( int i = 0; i < 3; i++ ) 
    {
      recoil->pos[i] = pka->pos[i]; 
      // used to assign the new position to the recoil, but
      // we have to make sure the recoil starts in the appropriate material!
      pka->pos[i] += pka->dir[i] * ( ls - simconf->tau );
      recoil->dir[i] = pka->dir[i] * p1;
    }
    recoil->t = pka->t;
    recoil->e = den;
    // displacement energy...
/*    if( recoil->z1 == 54 )  
      recoil->e -= 5.0; 
    else 
      recoil->e -= 25.0;*/
    recoil->m1 = material->element[nn]->m;
    recoil->z1 = material->element[nn]->z;

    // create a random vector perpendicular to pka.dir
    // there is a cleverer way by using the azimuthal angle of scatter...
    do
    { 
      for( int i = 0; i < 3; i++ ) rdir[i] = dr250() - 0.5;
      v_cross( pka->dir, rdir, perp );
      norm = sqrtf( v_dot( perp, perp) );
    }
    while( norm == 0.0 );
    v_scale( perp, 1.0 / norm );

    psi = atan( st / ( ct + material->element[nn]->my ) );
    v_scale( pka->dir, cos( psi ) );

    // calculate new direction, subtract from old dir (stored in recoil)
    for( int i = 0; i < 3; i++ ) 
    {
      pka->dir[i] += perp[i] * sin( psi );
      recoil->dir[i] -= pka->dir[i] * p2;
    }
    
    // end cascade if a CUT boundary is crossed
    terminate = false;
    for( int i = 0; i < 3; i++ ) 
      if( sample->bc[i] == sampleBase::CUT && ( pka->pos[i] > sample->w[i] || pka->pos[i] < 0.0 ) ) terminate = true;
    
    // if recoil energy > 100.0 eV, put the recoil on the stack
    if( spawnRecoil() && !terminate )
    {
      v_norm( recoil->dir );
      recoil->tag = material->tag;
      recoil->gen = pka->gen + 1;
      if( pka->md > 0 ) 
        recoil->md = pka->md +1;
      else 
        recoil->md = 0;
      recoil->id = simconf->id++;
      recoils.push( recoil );
      if( simconf->fullTraj ) printf( "spawn %d %d\n", recoil->id, pka->id );
    }
    else delete recoil;

    // output the full trajectory
    if( simconf->fullTraj )
      printf( "cont %f %f %f %d %d %d\n", pka->pos[0], pka->pos[1], pka->pos[2], pka->z1, pka->md, pka->id );

  } while ( pka->e > pka->ef && !terminate );

  if( simconf->fullTraj )
   if(pka->z1 == 54 && pka->gen > 0 ) printf( "\n" );
}


/*
materialBase* sampleType::lookupLayer( const float* pos ) 
{ 
  float dif[3];

  dif[0] = pos[0] - 100.0;
  dif[1] = pos[1];
  dif[2] = pos[2];
  float r2 = v_dot( dif, dif );
  if( r2 < 2500.0 ) // r<50.0
    return material[1]; 
  else
    return material[0]; 
}
*/
