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

/*! \file src/Healpix_2.15a/cxxsupport/pointing.h
 *  Class representing a direction in 3D space
 *
 *  Copyright (C) 2003 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_POINTING_H
#define PLANCK_POINTING_H

#include <cmath>
#include "src/Healpix_2.15a/vec3.h"
#include "src/Healpix_2.15a/cxxutils.h"

/*! \defgroup pointinggroup Pointings */
/*! \{ */

/*! Class representing a direction in 3D space or a location on the
    unit sphere. */
class pointing
  {
  public:
    /*! Colatitude of the pointing (i.e. the North pole is at \a theta=0). */
    double theta;
    /*! Longitude of the pointing. */
    double phi;

    /*! Default constructor. \a theta and \a phi are not initialized. */
    pointing() {}
    /*! Creates a pointing with \a Theta and \a Phi. */
    pointing (double Theta, double Phi) : theta(Theta), phi(Phi) {}

// FIXME: should become "explicit" some time
    /*! Creates a pointing from the vector \a inp. \a inp need not be
        normalized. */
    pointing (const vec3 &inp)
      {
      using namespace std;
      theta = atan2(sqrt(inp.x*inp.x+inp.y*inp.y),inp.z);
      phi = safe_atan2 (inp.y,inp.x);
      if (phi<0) phi += twopi;
      }
// FIXME: should be removed some time
    /*! Returns a normalized vector pointing in the same direction. */
    operator vec3() const
      {
      double st=sin(theta);
      return vec3 (st*cos(phi), st*sin(phi), cos(theta));
      }
    /*! Returns a normalized vector pointing in the same direction. */
    vec3 to_vec3() const
      {
      double st=sin(theta);
      return vec3 (st*cos(phi), st*sin(phi), cos(theta));
      }
    /*! Changes the angles so that \a 0<=theta<=pi and \a 0<=phi<2*pi. */
    void normalize()
      {
      theta=fmodulo(theta,twopi);
      if (theta>pi)
        {
        phi+=pi;
        theta=twopi-theta;
        }
      phi=fmodulo(phi,twopi);
      }
  };

/*! Converts \a vec to \a ptg. \a vec need not be normalized.
    \relates pointing */
inline void vec2pnt(const vec3 &vec, pointing &ptg)
  {
  using namespace std;
  ptg.theta = atan2(sqrt(vec.x*vec.x+vec.y*vec.y),vec.z);
  ptg.phi = safe_atan2 (vec.y,vec.x);
  if (ptg.phi<0) ptg.phi += twopi;
  }

/*! Writes \a p to \a os.
    \relates pointing */
inline std::ostream &operator<< (std::ostream &os, const pointing &p)
  {
  os << p.theta << ", " << p.phi << std::endl;
  return os;
  }

/*! \} */

#endif
