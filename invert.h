#ifndef INVERT_H
#define INVERT_H 1

#include <math.h>

class inverter
{
protected :
  float maxx, maxf, tol;

public :
  virtual float f( float x ) = 0;
  float x( float f );

  inverter() { maxx = 0.0; }
};


class massInverter : public inverter
{
public:
  virtual float f( float x );

  massInverter() { maxx = 235.0; tol = 1e-7; maxf = f( maxx ); }
};

class energyInverter : public inverter
{
 float A;
public:
  virtual float f( float x );

  energyInverter() { maxx = 186.98; tol = 1e-7; maxf = f( maxx ); }
};

#endif
