#ifndef INVERT_H
#define INVERT_H 1

#include <math.h>

class inverter
{
protected :
  double maxx, maxf, tol;

public :
  virtual double f( double x ) = 0;
  double x( double f );

  inverter() { maxx = 0.0; }
};


class massInverter : public inverter
{
public:
  virtual double f( double x );

  massInverter() { maxx = 235.0; tol = 1e-7; maxf = f( maxx ); }
};

class energyInverter : public inverter
{
 double A;
public:
  virtual double f( double x );

  energyInverter() { maxx = 186.98; tol = 1e-7; setMass(100.0); }
  void setMass( double _A ) {
    A = _A;
    maxf = f( maxx );
  }
};

#endif
