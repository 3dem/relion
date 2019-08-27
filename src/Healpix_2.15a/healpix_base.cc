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

/*
 *  Copyright (C) 2003, 2004, 2005, 2006 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "src/Healpix_2.15a/healpix_base.h"
#include "src/Healpix_2.15a/cxxutils.h"
#include "src/Healpix_2.15a/pointing.h"
#include "src/Healpix_2.15a/arr.h"
#include "src/Healpix_2.15a/geom_utils.h"

using namespace std;

short Healpix_Base::ctab[];
short Healpix_Base::utab[];

const nside_dummy SET_NSIDE=nside_dummy();

Healpix_Base::Tablefiller::Tablefiller()
  {
  for (int m=0; m<0x100; ++m)
    {
    ctab[m] =
         (m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
      | ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4);
    utab[m] =
         (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
      | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7);
    }
  }

Healpix_Base::Tablefiller Healpix_Base::Filler;

const int Healpix_Base::jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
const int Healpix_Base::jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

int Healpix_Base::npix2nside (int npix)
  {
  int res=isqrt(npix/12);
  planck_assert (npix==res*res*12, "npix2nside: invalid argument");
  return res;
  }

int Healpix_Base::ring_above (double z) const
  {
  double az=abs(z);
  if (az>twothird) // polar caps
    {
    int iring = int(nside_*sqrt(3*(1-az)));
    return (z>0) ? iring : 4*nside_-iring-1;
    }
  else // ----- equatorial region ---------
    return int(nside_*(2-1.5*z));
  }

void Healpix_Base::in_ring(int iz, double phi0, double dphi,
  vector<int> &listir) const
  {
  int nr, ir, ipix1;
  double shift=0.5;

  if (iz<nside_) // north pole
    {
    ir = iz;
    nr = ir*4;
    ipix1 = 2*ir*(ir-1);        //    lowest pixel number in the ring
    }
  else if (iz>(3*nside_)) // south pole
    {
    ir = 4*nside_ - iz;
    nr = ir*4;
    ipix1 = npix_ - 2*ir*(ir+1); // lowest pixel number in the ring
    }
  else // equatorial region
    {
    ir = iz - nside_ + 1;           //    within {1, 2*nside + 1}
    nr = nside_*4;
    if ((ir&1)==0) shift = 0;
    ipix1 = ncap_ + (ir-1)*nr; // lowest pixel number in the ring
    }

  int ipix2 = ipix1 + nr - 1;       //    highest pixel number in the ring

   // ----------- constructs the pixel list --------------
  if (dphi > (pi-1e-7))
    for (int i=ipix1; i<=ipix2; ++i) listir.push_back(i);
  else
    {
    int ip_lo = ifloor<int>(nr*inv_twopi*(phi0-dphi) - shift)+1;
    int ip_hi = ifloor<int>(nr*inv_twopi*(phi0+dphi) - shift);
    int pixnum = ip_lo+ipix1;
    if (pixnum<ipix1) pixnum += nr;
    for (int i=ip_lo; i<=ip_hi; ++i, ++pixnum)
      {
      if (pixnum>ipix2) pixnum -= nr;
      listir.push_back(pixnum);
      }
    }
  }

void Healpix_Base::nest2xyf (int pix, int &ix, int &iy, int &face_num) const
  {
  face_num = pix>>(2*order_);
  pix &= (npface_-1);
  int raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  ix = ctab[raw&0xff] | (ctab[raw>>8]<<4);
  pix >>= 1;
  raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  iy = ctab[raw&0xff] | (ctab[raw>>8]<<4);
  }

int Healpix_Base::xyf2nest (int ix, int iy, int face_num) const
  {
  return (face_num<<(2*order_)) +
      (utab[ix&0xff] | (utab[ix>>8]<<16)
    | (utab[iy&0xff]<<1) | (utab[iy>>8]<<17));
  }

void Healpix_Base::ring2xyf (int pix, int &ix, int &iy, int &face_num) const
  {
  int iring, iphi, kshift, nr;

  int nl2 = 2*nside_;

  if (pix<ncap_) // North Polar cap
    {
    iring = int(0.5*(1+isqrt(1+2*pix))); //counted from North pole
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    face_num=0;
    int tmp = iphi-1;
    if (tmp>=(2*iring))
      {
      face_num=2;
      tmp-=2*iring;
      }
    if (tmp>=iring) ++face_num;
    }
  else if (pix<(npix_-ncap_)) // Equatorial region
    {
    int ip = pix - ncap_;
    if (order_>=0)
      {
      iring = (ip>>(order_+2)) + nside_; // counted from North pole
      iphi  = (ip&(4*nside_-1)) + 1;
      }
    else
      {
      iring = (ip/(4*nside_)) + nside_; // counted from North pole
      iphi  = (ip%(4*nside_)) + 1;
      }
    kshift = (iring+nside_)&1;
    nr = nside_;
    unsigned int ire = iring-nside_+1;
    unsigned int irm = nl2+2-ire;
    int ifm, ifp;
    if (order_>=0)
      {
      ifm = (iphi - ire/2 + nside_ -1) >> order_;
      ifp = (iphi - irm/2 + nside_ -1) >> order_;
      }
    else
      {
      ifm = (iphi - ire/2 + nside_ -1) / nside_;
      ifp = (iphi - irm/2 + nside_ -1) / nside_;
      }
    if (ifp == ifm) // faces 4 to 7
      face_num = (ifp==4) ? 4 : ifp+4;
    else if (ifp<ifm) // (half-)faces 0 to 3
      face_num = ifp;
    else // (half-)faces 8 to 11
      face_num = ifm + 8;
    }
  else // South Polar cap
    {
    int ip = npix_ - pix;
    iring = int(0.5*(1+isqrt(2*ip-1))); //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    face_num=8;
    int tmp = iphi-1;
    if (tmp>=(2*nr))
      {
      face_num=10;
      tmp-=2*nr;
      }
    if (tmp>=nr) ++face_num;
    }

  int irt = iring - (jrll[face_num]*nside_) + 1;
  int ipt = 2*iphi- jpll[face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  ix =  (ipt-irt) >>1;
  iy =(-(ipt+irt))>>1;
  }

int Healpix_Base::xyf2ring (int ix, int iy, int face_num) const
  {
  int nl4 = 4*nside_;
  int jr = (jrll[face_num]*nside_) - ix - iy  - 1;

  int nr, kshift, n_before;
  if (jr<nside_)
    {
    nr = jr;
    n_before = 2*nr*(nr-1);
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    n_before = npix_ - 2*(nr+1)*nr;
    kshift = 0;
    }
  else
    {
    nr = nside_;
    n_before = ncap_ + (jr-nside_)*nl4;
    kshift = (jr-nside_)&1;
    }

  int jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;

  return n_before + jp - 1;
  }

double Healpix_Base::ring2z (int ring) const
  {
  if (ring<nside_)
    return 1 - ring*ring*fact2_;
  if (ring <=3*nside_)
    return (2*nside_-ring)*fact1_;
  ring=4*nside_ - ring;
  return ring*ring*fact2_ - 1;
  }

int Healpix_Base::pix2ring (int pix) const
  {
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      return int(0.5*(1+isqrt(1+2*pix))); //counted from North pole
    else if (pix<(npix_-ncap_)) // Equatorial region
      {
      int ip  = pix - ncap_;
      return ip/(4*nside_) + nside_; // counted from North pole
      }
    else // South Polar cap
      {
      int ip = npix_ - pix;
      return 4*nside_ - int(0.5*(1+isqrt(2*ip-1))); //counted from South pole
      }
    }
  else
    {
    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);
    return (jrll[face_num]<<order_) - ix - iy - 1;
    }
  }

int Healpix_Base::nest2ring (int pix) const
  {
  planck_assert(order_>=0, "nest2ring: need hierarchical map");
  int ix, iy, face_num;
  nest2xyf (pix, ix, iy, face_num);
  return xyf2ring (ix, iy, face_num);
  }

int Healpix_Base::ring2nest (int pix) const
  {
  planck_assert(order_>=0, "ring2nest: need hierarchical map");
  int ix, iy, face_num;
  ring2xyf (pix, ix, iy, face_num);
  return xyf2nest (ix, iy, face_num);
  }

int Healpix_Base::nest2peano (int pix) const
  {
  static const unsigned char subpix[8][4] = {
    { 0, 1, 3, 2 }, { 3, 0, 2, 1 }, { 2, 3, 1, 0 }, { 1, 2, 0, 3 },
    { 0, 3, 1, 2 }, { 1, 0, 2, 3 }, { 2, 1, 3, 0 }, { 3, 2, 0, 1 } };
  const unsigned char subpath[8][4] = {
    { 4, 0, 6, 0 }, { 7, 5, 1, 1 }, { 2, 4, 2, 6 }, { 3, 3, 7, 5 },
    { 0, 2, 4, 4 }, { 5, 1, 5, 3 }, { 6, 6, 0, 2 }, { 1, 7, 3, 7 } };
  static const unsigned char face2path[12] = {
    2, 5, 2, 5, 3, 6, 3, 6, 2, 3, 2, 3 };
  static const unsigned char face2peanoface[12] = {
    0, 5, 6, 11, 10, 1, 4, 7, 2, 3, 8, 9 };

  int face = pix>>(2*order_);
  unsigned char path = face2path[face];
  int result = 0;

  for (int shift=2*order_-2; shift>=0; shift-=2)
    {
    unsigned char spix = (pix>>shift) & 0x3;
    result <<= 2;
    result |= subpix[path][spix];
    path=subpath[path][spix];
    }

  return result + (int(face2peanoface[face])<<(2*order_));
  }

int Healpix_Base::peano2nest (int pix) const
  {
  static const unsigned char subpix[8][4] = {
    { 0, 1, 3, 2 }, { 1, 3, 2, 0 }, { 3, 2, 0, 1 }, { 2, 0, 1, 3 },
    { 0, 2, 3, 1 }, { 1, 0, 2, 3 }, { 3, 1, 0, 2 }, { 2, 3, 1, 0 } };
  static const unsigned char subpath[8][4] = {
    { 4, 0, 0, 6 }, { 5, 1, 1, 7 }, { 6, 2, 2, 4 }, { 7, 3, 3, 5 },
    { 0, 4, 4, 2 }, { 1, 5, 5, 3 }, { 2, 6, 6, 0 }, { 3, 7, 7, 1 } };
  static const unsigned char face2path[12] = {
    2, 6, 2, 3, 3, 5, 2, 6, 2, 3, 3, 5 };
  static const unsigned char peanoface2face[12] = {
    0, 5, 8, 9, 6, 1, 2, 7, 10, 11, 4, 3 };

  int face = pix>>(2*order_);
  unsigned char path = face2path[face];
  int result = 0;

  for (int shift=2*order_-2; shift>=0; shift-=2)
    {
    unsigned char spix = (pix>>shift) & 0x3;
    result <<= 2;
    result |= subpix[path][spix];
    path=subpath[path][spix];
    }

  return result + (int(peanoface2face[face])<<(2*order_));
  }

int Healpix_Base::ang2pix_z_phi (double z, double phi) const
  {
  double za = abs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; // in [0,4)

  if (scheme_==RING)
    {
    if (za<=twothird) // Equatorial region
      {
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*z*0.75;
      int jp = int(temp1-temp2); // index of  ascending edge line
      int jm = int(temp1+temp2); // index of descending edge line

      // ring number counted from z=2/3
      int ir = nside_ + 1 + jp - jm; // in {1,2n+1}
      int kshift = 1-(ir&1); // kshift=1 if ir even, 0 otherwise

      int ip = (jp+jm-nside_+kshift+1)/2; // in {0,4n-1}
      ip = imodulo(ip,4*nside_);

      return ncap_ + (ir-1)*4*nside_ + ip;
      }
    else  // North & South polar caps
      {
      double tp = tt-int(tt);
      double tmp = nside_*sqrt(3*(1-za));

      int jp = int(tp*tmp); // increasing edge line index
      int jm = int((1.0-tp)*tmp); // decreasing edge line index

      int ir = jp+jm+1; // ring number counted from the closest pole
      int ip = int(tt*ir); // in {0,4*ir-1}
      ip = imodulo(ip,4*ir);

      if (z>0)
        return 2*ir*(ir-1) + ip;
      else
        return npix_ - 2*ir*(ir+1) + ip;
      }
    }
  else // scheme_ == NEST
    {
    int face_num, ix, iy;

    if (za<=twothird) // Equatorial region
      {
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*(z*0.75);
      int jp = int(temp1-temp2); // index of  ascending edge line
      int jm = int(temp1+temp2); // index of descending edge line
      int ifp = jp >> order_;  // in {0,4}
      int ifm = jm >> order_;
      if (ifp == ifm)           // faces 4 to 7
        face_num = (ifp==4) ? 4: ifp+4;
      else if (ifp < ifm)       // (half-)faces 0 to 3
        face_num = ifp;
      else                      // (half-)faces 8 to 11
        face_num = ifm + 8;

      ix = jm & (nside_-1);
      iy = nside_ - (jp & (nside_-1)) - 1;
      }
    else // polar region, za > 2/3
      {
      int ntt = int(tt);
      if (ntt>=4) ntt=3;
      double tp = tt-ntt;
      double tmp = nside_*sqrt(3*(1-za));

      int jp = int(tp*tmp); // increasing edge line index
      int jm = int((1.0-tp)*tmp); // decreasing edge line index
      if (jp>=nside_) jp = nside_-1; // for points too close to the boundary
      if (jm>=nside_) jm = nside_-1;
      if (z >= 0)
        {
        face_num = ntt;  // in {0,3}
        ix = nside_ - jm - 1;
        iy = nside_ - jp - 1;
        }
      else
        {
        face_num = ntt + 8; // in {8,11}
        ix =  jp;
        iy =  jm;
        }
      }

    return xyf2nest(ix,iy,face_num);
    }
  }

void Healpix_Base::pix2ang_z_phi (int pix, double &z, double &phi) const
  {
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      {
      int iring = int(0.5*(1+isqrt(1+2*pix))); //counted from North pole
      int iphi  = (pix+1) - 2*iring*(iring-1);

      z = 1.0 - (iring*iring)*fact2_;
      phi = (iphi-0.5) * halfpi/iring;
      }
    else if (pix<(npix_-ncap_)) // Equatorial region
      {
      int ip  = pix - ncap_;
      int iring = ip/(4*nside_) + nside_; // counted from North pole
      int iphi  = ip%(4*nside_) + 1;
      // 1 if iring+nside is odd, 1/2 otherwise
      double fodd = ((iring+nside_)&1) ? 1 : 0.5;

      int nl2 = 2*nside_;
      z = (nl2-iring)*fact1_;
      phi = (iphi-fodd) * pi/nl2;
      }
    else // South Polar cap
      {
      int ip = npix_ - pix;
      int iring = int(0.5*(1+isqrt(2*ip-1))); //counted from South pole
      int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

      z = -1.0 + (iring*iring)*fact2_;
      phi = (iphi-0.5) * halfpi/iring;
      }
    }
  else
    {
    int nl4 = nside_*4;

    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);

    int jr = (jrll[face_num]<<order_) - ix - iy - 1;

    int nr, kshift;
    if (jr<nside_)
      {
      nr = jr;
      z = 1 - nr*nr*fact2_;
      kshift = 0;
      }
    else if (jr > 3*nside_)
      {
      nr = nl4-jr;
      z = nr*nr*fact2_ - 1;
      kshift = 0;
      }
    else
      {
      nr = nside_;
      z = (2*nside_-jr)*fact1_;
      kshift = (jr-nside_)&1;
      }

    int jp = (jpll[face_num]*nr + ix -iy + 1 + kshift) / 2;
    if (jp>nl4) jp-=nl4;
    if (jp<1) jp+=nl4;

    phi = (jp-(kshift+1)*0.5)*(halfpi/nr);
    }
  }

void Healpix_Base::query_disc (const pointing &ptg, double radius,
  vector<int>& listpix) const
  {
  listpix.clear();

  double dth1 = fact2_;
  double dth2 = fact1_;
  double cosang = cos(radius);

  double z0 = cos(ptg.theta);
  double xa = 1./sqrt((1-z0)*(1+z0));

  double rlat1  = ptg.theta - radius;
  double zmax = cos(rlat1);
  int irmin = ring_above (zmax)+1;

  if (rlat1<=0) // north pole in the disc
    for (int m=1; m<irmin; ++m) // rings completely in the disc
      in_ring (m, 0, pi, listpix);

  double rlat2  = ptg.theta + radius;
  double zmin = cos(rlat2);
  int irmax = ring_above (zmin);

// ------------- loop on ring number ---------------------
  for (int iz=irmin; iz<=irmax; ++iz) // rings partially in the disc
    {
    double z;
    if (iz<nside_) // north polar cap
      z = 1.0 - iz*iz*dth1;
    else if (iz <= (3*nside_)) // tropical band + equat.
      z = (2*nside_-iz) * dth2;
    else
      z = -1.0 + (4*nside_-iz)*(4*nside_-iz)*dth1;

// --------- phi range in the disc for each z ---------
    double x = (cosang-z*z0)*xa;
    double ysq = 1-z*z-x*x;
    planck_assert(ysq>=0, "error in query_disc()");
    double dphi=atan2(sqrt(ysq),x);
    in_ring (iz, ptg.phi, dphi, listpix);
    }

  if (rlat2>=pi) // south pole in the disc
    for (int m=irmax+1; m<(4*nside_); ++m)  // rings completely in the disc
      in_ring (m, 0, pi, listpix);

  if (scheme_==NEST)
    for (unsigned int m=0; m<listpix.size(); ++m)
      listpix[m] = ring2nest(listpix[m]);
  }

void Healpix_Base::get_ring_info (int ring, int &startpix, int &ringpix,
  double &costheta, double &sintheta, bool &shifted) const
  {
  planck_assert(scheme_==RING,"map must be in RING scheme");
  int northring = (ring>2*nside_) ? 4*nside_-ring : ring;
  if (northring < nside_)
    {
    double tmp = northring*northring*fact2_;
    costheta = 1 - tmp;
    sintheta = sqrt(tmp*(2-tmp));
    ringpix = 4*northring;
    shifted = true;
    startpix = 2*northring*(northring-1);
    }
  else
    {
    costheta = (2*nside_-northring)*fact1_;
    sintheta = sqrt((1+costheta)*(1-costheta));
    ringpix = 4*nside_;
    shifted = ((northring-nside_) & 1) == 0;
    startpix = ncap_ + (northring-nside_)*ringpix;
    }
  if (northring != ring) // southern hemisphere
    {
    costheta = -costheta;
    startpix = npix_ - startpix - ringpix;
    }
  }

void Healpix_Base::neighbors (int pix, fix_arr<int,8> &result) const
  {
  static const int xoffset[] = { -1,-1, 0, 1, 1, 1, 0,-1 };
  static const int yoffset[] = {  0, 1, 1, 1, 0,-1,-1,-1 };
  static const int facearray[][12] =
        { {  8, 9,10,11,-1,-1,-1,-1,10,11, 8, 9 },   // S
          {  5, 6, 7, 4, 8, 9,10,11, 9,10,11, 8 },   // SE
          { -1,-1,-1,-1, 5, 6, 7, 4,-1,-1,-1,-1 },   // E
          {  4, 5, 6, 7,11, 8, 9,10,11, 8, 9,10 },   // SW
          {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 },   // center
          {  1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4 },   // NE
          { -1,-1,-1,-1, 7, 4, 5, 6,-1,-1,-1,-1 },   // W
          {  3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7 },   // NW
          {  2, 3, 0, 1,-1,-1,-1,-1, 0, 1, 2, 3 } }; // N
  static const int swaparray[][12] =
        { {  0,0,0,0,0,0,0,0,3,3,3,3 },   // S
          {  0,0,0,0,0,0,0,0,6,6,6,6 },   // SE
          {  0,0,0,0,0,0,0,0,0,0,0,0 },   // E
          {  0,0,0,0,0,0,0,0,5,5,5,5 },   // SW
          {  0,0,0,0,0,0,0,0,0,0,0,0 },   // center
          {  5,5,5,5,0,0,0,0,0,0,0,0 },   // NE
          {  0,0,0,0,0,0,0,0,0,0,0,0 },   // W
          {  6,6,6,6,0,0,0,0,0,0,0,0 },   // NW
          {  3,3,3,3,0,0,0,0,0,0,0,0 } }; // N

  int ix, iy, face_num;
  (scheme_==RING) ?
    ring2xyf(pix,ix,iy,face_num) : nest2xyf(pix,ix,iy,face_num);

  const int nsm1 = nside_-1;
  if ((ix>0)&&(ix<nsm1)&&(iy>0)&&(iy<nsm1))
    {
    if (scheme_==RING)
      for (int m=0; m<8; ++m)
        result[m] = xyf2ring(ix+xoffset[m],iy+yoffset[m],face_num);
    else
      for (int m=0; m<8; ++m)
        result[m] = xyf2nest(ix+xoffset[m],iy+yoffset[m],face_num);
    }
  else
    {
    for (int i=0; i<8; ++i)
      {
      int x=ix+xoffset[i];
      int y=iy+yoffset[i];
      int nbnum=4;
      if (x<0)
        { x+=nside_; nbnum-=1; }
      else if (x>=nside_)
        { x-=nside_; nbnum+=1; }
      if (y<0)
        { y+=nside_; nbnum-=3; }
      else if (y>=nside_)
        { y-=nside_; nbnum+=3; }

      int f = facearray[nbnum][face_num];
      if (f>=0)
        {
        if (swaparray[nbnum][face_num]&1) x=nside_-x-1;
        if (swaparray[nbnum][face_num]&2) y=nside_-y-1;
        if (swaparray[nbnum][face_num]&4) std::swap(x,y);
        result[i] = (scheme_==RING) ? xyf2ring(x,y,f) : xyf2nest(x,y,f);
        }
      else
        result[i] = -1;
      }
    }
  }

void Healpix_Base::get_ring_info2 (int ring, int &startpix, int &ringpix,
  double &theta, bool &shifted) const
  {
  int northring = (ring>2*nside_) ? 4*nside_-ring : ring;
  if (northring < nside_)
    {
    double tmp = northring*northring*fact2_;
    double costheta = 1 - tmp;
    double sintheta = sqrt(tmp*(2-tmp));
    theta = atan2(sintheta,costheta);
    ringpix = 4*northring;
    shifted = true;
    startpix = 2*northring*(northring-1);
    }
  else
    {
    theta = acos((2*nside_-northring)*fact1_);
    ringpix = 4*nside_;
    shifted = ((northring-nside_) & 1) == 0;
    startpix = ncap_ + (northring-nside_)*ringpix;
    }
  if (northring != ring) // southern hemisphere
    {
    theta = pi-theta;
    startpix = npix_ - startpix - ringpix;
    }
  }

void Healpix_Base::get_interpol (const pointing &ptg, fix_arr<int,4> &pix,
  fix_arr<double,4> &wgt) const
  {
  double z = cos (ptg.theta);
  int ir1 = ring_above(z);
  int ir2 = ir1+1;
  double theta1, theta2, w1, tmp, dphi;
  int sp,nr;
  bool shift;
  int i1,i2;
  if (ir1>0)
    {
    get_ring_info2 (ir1, sp, nr, theta1, shift);
    dphi = twopi/nr;
    tmp = (ptg.phi/dphi - .5*shift);
    i1 = (tmp<0) ? int(tmp)-1 : int(tmp);
    w1 = (ptg.phi-(i1+.5*shift)*dphi)/dphi;
    i2 = i1+1;
    if (i1<0) i1 +=nr;
    if (i2>=nr) i2 -=nr;
    pix[0] = sp+i1; pix[1] = sp+i2;
    wgt[0] = 1-w1; wgt[1] = w1;
    }
  if (ir2<(4*nside_))
    {
    get_ring_info2 (ir2, sp, nr, theta2, shift);
    dphi = twopi/nr;
    tmp = (ptg.phi/dphi - .5*shift);
    i1 = (tmp<0) ? int(tmp)-1 : int(tmp);
    w1 = (ptg.phi-(i1+.5*shift)*dphi)/dphi;
    i2 = i1+1;
    if (i1<0) i1 +=nr;
    if (i2>=nr) i2 -=nr;
    pix[2] = sp+i1; pix[3] = sp+i2;
    wgt[2] = 1-w1; wgt[3] = w1;
    }

  if (ir1==0)
    {
    double wtheta = ptg.theta/theta2;
    wgt[2] *= wtheta; wgt[3] *= wtheta;
    double fac = (1-wtheta)*0.25;
    wgt[0] = fac; wgt[1] = fac; wgt[2] += fac; wgt[3] +=fac;
    pix[0] = (pix[2]+2)%4;
    pix[1] = (pix[3]+2)%4;
    }
  else if (ir2==4*nside_)
    {
    double wtheta = (ptg.theta-theta1)/(pi-theta1);
    wgt[0] *= (1-wtheta); wgt[1] *= (1-wtheta);
    double fac = wtheta*0.25;
    wgt[0] += fac; wgt[1] += fac; wgt[2] = fac; wgt[3] =fac;
    pix[2] = ((pix[0]+2)&3)+npix_-4;
    pix[3] = ((pix[1]+2)&3)+npix_-4;
    }
  else
    {
    double wtheta = (ptg.theta-theta1)/(theta2-theta1);
    wgt[0] *= (1-wtheta); wgt[1] *= (1-wtheta);
    wgt[2] *= wtheta; wgt[3] *= wtheta;
    }

  if (scheme_==NEST)
    for (int m=0; m<pix.size(); ++m)
      pix[m] = ring2nest(pix[m]);
  }

int Healpix_Base::get_npix() const
{
	return npix_;
}

double Healpix_Base::max_pixrad() const
  {
  vec3 va,vb;
  va.set_z_phi (2./3., pi/(4*nside_));
  double t1 = 1.-1./nside_;
  t1*=t1;
  vb.set_z_phi (1-t1/3, 0);
  return v_angle(va,vb);
  }
