#ifndef MATERIAL_H
#define MATERIAL_H 1

#include <vector>
#include <stdio.h>

#include "ion.h"
#include "element.h"

void v_cross( const float *a1, const float *a2, float *b );
void v_scale( float *a1, float b ); // in=place scale
float v_dot( const float *a1, const float *a2 );
void v_norm( float *a1, float b = 1.0 ); // in-place normalize to 1 (or b)

inline float sqr( float a ) { return a*a; }
inline float cub( float a ) { return a*a*a; }

using namespace std;

struct materialBase {
  float rho;

  float am, az; // average mass and atomic number
  float arho, mu;
  float a, f, epsdg;
  float fd, kd, pmax;

  int tag;

  vector<elementBase*> element;

  //layerType() { semax = 0.0; sem = 0.0; sez = 0; }
  materialBase( float _rho ) : rho(_rho) {};

  // make sure stoiciometry is normalized, compute averages independent of pka
  void prepare();

  // compute pka dependent averages 
  void average( const ionBase *pka );
  float getrstop( const ionBase *pka );

protected:
  float rpstop( int z2, float e );
  float rstop( const ionBase *ion, int z2 );
};

#endif
