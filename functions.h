/*
MyTRIM - a three dimensional binary collision Monte Carlo library.
Copyright (C) 2008-2018  Daniel Schwen <daniel@schwen.de>

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA
*/

#ifndef MYTRIM_FUNCTIONS_H
#define MYTRIM_FUNCTIONS_H

#include <cmath>
#include <stdlib.h>

namespace MyTRIM_NS
{

inline void
v_cross(const Real * a1, const Real * a2, Real * b)
{
  for (unsigned int i = 0; i < 3; ++i)
    b[i] = a1[(i + 1) % 3] * a2[(i + 2) % 3] - a1[(i + 2) % 3] * a2[(i + 1) % 3];
}

inline void
v_cross(const Point & a1, const Point & a2, Point & b)
{
  for (unsigned int i = 0; i < 3; ++i)
    b(i) = a1((i + 1) % 3) * a2((i + 2) % 3) - a1((i + 2) % 3) * a2((i + 1) % 3);
}

inline void
v_scale(Real * a1, Real b) // in=place scale
{
  for (unsigned int i = 0; i < 3; ++i)
    a1[i] = a1[i] * b;
}

inline Real
v_dot(const Real * a1, const Real * a2)
{
  Real b = 0.0;
  for (int i = 0; i < 3; ++i)
    b += a1[i] * a2[i];
  return b;
}

inline void
v_norm(Real * a1, Real b = 1.0) // in-place normalize to b (= 1.0 default)
{
  v_scale(a1, b / std::sqrt(v_dot(a1, a1)));
}

inline void
v_norm(Point & a1, Real b = 1.0) // in-place normalize to b (= 1.0 default)
{
  a1 *= b / a1.norm();
}

inline Real
sqr(Real a)
{
  return a * a;
}
inline Real
cub(Real a)
{
  return a * a * a;
}
} // namespace MyTRIM_NS

#endif
