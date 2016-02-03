/*
MyTRIM - a three dimensional binary collision Monte Carlo library.
Copyright (C) 2008-2015  Daniel Schwen <daniel@schwen.de>

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

#include "point.h"

Point::Point()
{
  data[0] = 0.0;
  data[1] = 0.0;
  data[2] = 0.0;
}

Point::Point(Real x, Real y, Real z)
{
  data[0] = x;
  data[1] = y;
  data[2] = z;
}

Real &
Point::operator() (unsigned int i)
{
  return data[i];
}

const Real &
Point::operator() (unsigned int i) const
{
  return data[i];
}


Point
Point::operator+ (const Point & rhs)
{
  return Point(data[0] + rhs.data[0],
               data[1] + rhs.data[1],
               data[2] + rhs.data[2]);
}

Point
Point::operator- (const Point & rhs)
{
  return Point(data[0] - rhs.data[0],
               data[1] - rhs.data[1],
               data[2] - rhs.data[2]);
}

Point
Point::operator* (Real rhs)
{
  return Point(data[0] * rhs,
               data[1] * rhs,
               data[2] * rhs);
}

Point
Point::operator/ (Real rhs)
{
  return Point(data[0] / rhs,
               data[1] / rhs,
               data[2] / rhs);
}


Point
Point::operator- ()
{
  return Point(-data[0],
               -data[1],
               -data[2]);
}


Point &
Point::operator+=(const Point & rhs)
{
  data[0] += rhs.data[0];
  data[1] += rhs.data[1];
  data[2] += rhs.data[2];
  return *this;
}

Point &
Point::operator-=(const Point & rhs)
{
  data[0] -= rhs.data[0];
  data[1] -= rhs.data[1];
  data[2] -= rhs.data[2];
  return *this;
}

Point &
Point::operator*=(Real rhs)
{
  data[0] *= rhs;
  data[1] *= rhs;
  data[2] *= rhs;
  return *this;
}

Point &
Point::operator/=(Real rhs)
{
  data[0] /= rhs;
  data[1] /= rhs;
  data[2] /= rhs;
  return *this;
}
