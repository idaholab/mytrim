#ifndef MYTRIM_FUNCTIONS_H
#define MYTRIM_FUNCTIONS_H

#include <cmath>
#include <stdlib.h>

namespace MyTRIM_NS {

inline void v_cross(const Real * a1, const Real * a2, Real * b)
{
  for (unsigned int i = 0; i < 3; ++i)
    b[i] = a1[(i+1)%3] * a2[(i+2)%3] - a1[(i+2)%3] * a2[(i+1)%3];
}

inline void v_cross(const Point & a1, const Point & a2, Point & b)
{
  for (unsigned int i = 0; i < 3; ++i)
    b[i] = a1[(i+1)%3] * a2[(i+2)%3] - a1[(i+2)%3] * a2[(i+1)%3];
}

inline void v_scale(Real * a1, Real b) // in=place scale
{
  for (unsigned int i = 0; i < 3; ++i)
    a1[i] = a1[i] * b;
}

inline Real v_dot(const Real * a1, const Real * a2)
{
  Real b = 0.0;
  for (int i = 0; i < 3; i++) b += a1[i] * a2[i];
    return b;
}

inline void v_norm(Real * a1, Real b = 1.0) // in-place normalize to b (= 1.0 default)
{
  v_scale(a1, b / std::sqrt(v_dot(a1, a1)));
}

inline void v_norm(Point & a1, Real b = 1.0) // in-place normalize to b (= 1.0 default)
{
  a1 *= b / a1.size();
}

inline Real sqr(Real a) { return a*a; }
inline Real cub(Real a) { return a*a*a; }

}

#endif
