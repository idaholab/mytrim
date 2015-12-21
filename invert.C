#include "invert.h"

#include <stdio.h>

using namespace MyTRIM_NS;

Real inverter::x(Real f1)
{
  Real f2, x1 = maxx / 2.0;
  Real w = maxx / 4.0;

  // no point in doing more than 32 iterations for Real (precission)
  for (int i = 0; i < 32; i++)
  {
    f2 = f(x1) / maxf;
    if (std::abs(f2 - f1) <= tol) break;
    if (f2 > f1)
      x1 -= w;
    else
      x1 += w;
    w *= 0.5;
  }

  return x1;
}

Real massInverter::f(Real x)
{
  return (  100.088
           + 0.112798 * std::erf(-5.56257 + 0.0471405 * x)
           + 37.4781 * std::erf(-19.3772 + 0.137386 * x)
           + 37.4781 * std::erf(-13.0462 + 0.137386 * x)
           + 12.5094 * std::erf(-30.8853 + 0.229537 * x)
           + 12.5094 * std::erf(-23.2853 + 0.229537 * x)
        ) / 200.1756;
}

Real energyInverter::f(Real x)
{
  Real x1 = x / (1.0 - A / 234.0);
  return (-0.00014122 + (0.00014122 -7.12299E-7 * x1) * std::exp(0.0886603 * x1)) / 127.216;
}
