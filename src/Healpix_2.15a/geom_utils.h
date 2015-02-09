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

/*! \file src/Healpix_2.15a/cxxsupport/geom_utils.h
 *  Geometric utility functions.
 *
 *  Copyright (C) 2003, 2006 Max-Planck-Society
 *  \author Martin Reinecke
 *  \author Reinhard Hell
 */

#include "src/Healpix_2.15a/cxxutils.h"
#include "src/Healpix_2.15a/vec3.h"

/*! Returns the orientation when looking from point \a loc on the unit
    sphere in the direction \a dir. \a loc must be normalized. The result
    ranges from -pi to pi, is 0 for North and pi/2 for West, i.e. the angle
    is given in mathematically positive sense.

    If \a loc is the North or South pole, the returned angle is
    \a atan2(dir.y,dir.x). */
inline double orientation (const vec3 &loc, const vec3 &dir)
  {
// FIXME: here is still optimization potential
  if (loc.x==0 && loc.y==0)
    {
    if (loc.z>0) return safe_atan2(dir.y,-dir.x);
    else return safe_atan2(dir.y,dir.x);
    }
  vec3 east (-loc.y, loc.x, 0);
  vec3 north = crossprod(loc,east);
  double y = dotprod(dir,east);
  double x = dotprod(dir,north);
  return safe_atan2(-y,x);
  }

/*! Returns the angle between \a v1 and \a v2 in radians. */
inline double v_angle (const vec3 &v1, const vec3 &v2)
  {
  return atan2 (crossprod(v1,v2).Length(), dotprod(v1,v2));
  }
