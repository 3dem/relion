/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file src/Healpix_2.15a/cxxsupport/vec3.h
 *  Class representing 3D cartesian vectors
 *
 *  Copyright (C) 2003, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_VEC3_H
#define PLANCK_VEC3_H

#include <cmath>
#include <iostream>

/*! \defgroup vec3group 3D vectors */
/*! \{ */

/*! Class representing a 3D cartesian vector. */
class vec3
  {
  public:
    double x, /*!< x-coordinate */
           y, /*!< y-coordinate */
           z; /*!< z-coordinate */

    /*! Default constructor. Does not initialize \a x, \a y, and \a z. */
    vec3 () {}
    /*! Creates a vector with the coordinates \a xc, \a yc, and \a zc. */
    vec3 (double xc, double yc, double zc)
      : x(xc), y(yc), z(zc) {}

    /*! Creates a unit vector from a z coordinate and an azimuthal angle. */
    void set_z_phi (double z_, double phi)
      {
      using namespace std;
      double sintheta = sqrt((1.-z_)*(1.+z_));
      x = sintheta*cos(phi);
      y = sintheta*sin(phi);
      z = z_;
      }

    /*! Normalizes the vector to length 1. */
    void Normalize ()
      {
      using namespace std;
      double l = 1.0/sqrt (x*x + y*y + z*z);
      x*=l; y*=l; z*=l;
      }

    /*! Returns the length of the vector. */
    double Length () const
      { return sqrt (x*x + y*y + z*z); }

    /*! Returns the squared length of the vector. */
    double SquaredLength () const
      { return (x*x + y*y + z*z); }
    /*! Returns the vector with the signs of all coordinates flipped. */
    const vec3 operator- () const
      { return vec3 (-x, -y, -z); }
    /*! Flips the signs of all coordinates. */
    void Flip ()
      { x=-x; y=-y; z=-z; }
    /*! Subtracts \a vec from the vector. */
    const vec3 operator- (const vec3 &vec) const
      { return vec3 (x-vec.x, y-vec.y, z-vec.z); }
    /*! Adds \a vec to the vector. */
    const vec3 operator+ (const vec3 &vec) const
      { return vec3 (x+vec.x, y+vec.y, z+vec.z); }
    /*! Returns the vector scaled by \a fact. */
    const vec3 operator* (double fact) const
      { return vec3 (x*fact, y*fact, z*fact); }
    /*! Returns the vector scaled by \a 1/fact. */
    const vec3 operator/ (double fact) const
      { double xfact = 1./fact; return vec3 (x*xfact, y*xfact, z*xfact); }
    /*! Scales the vector by \a fact. */
    vec3 &operator*= (double fact)
      { x*=fact; y*=fact; z*=fact; return *this; }
  };

/*! Returns the dot product of \a v1 and \a v2.
    \relates vec3 */
inline double dotprod(const vec3 &v1, const vec3 &v2)
  { return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z; }

/*! Returns the cross product of \a a and \a b.
    \relates vec3 */
inline vec3 crossprod(const vec3 &a, const vec3 &b)
  { return vec3 (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }

/*! Writes \a v to \a os.
    \relates vec3 */
inline std::ostream &operator<< (std::ostream &os, const vec3 &v)
  {
  os << v.x << ", " << v.y << ", " << v.z << std::endl;
  return os;
  }

/*! \} */

#endif
