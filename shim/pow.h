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

This file contains GPL code from libMesh:
Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
*/

#ifndef MYTRIM_POW_H
#define MYTRIM_POW_H

namespace Utility
{

template <int N, typename T>
struct do_pow
{
  static inline T apply(const T & x)
  {
    if (N % 2) // odd exponent
      return x * do_pow<N - 1, T>::apply(x);

    const T xNover2 = do_pow<N / 2, T>::apply(x);

    return xNover2 * xNover2;
  }
};

template <typename T>
struct do_pow<6, T>
{
  static inline T apply(const T & x)
  {
    const T x2 = x * x, x4 = x2 * x2;

    return x4 * x2;
  }
};

template <typename T>
struct do_pow<1, T>
{
  static inline T apply(const T & x) { return x; }
};

template <typename T>
struct do_pow<0, T>
{
  static inline T apply(const T &) { return 1; }
};

template <int N, typename T>
inline T
pow(const T & x)
{
  return do_pow<N, T>::apply(x);
}

} // namespace Utility

#endif // MYTRIM_POW_H
