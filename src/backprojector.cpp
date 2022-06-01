/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/*
 * backprojector.cpp
 *
 *  Created on: 24 Aug 2010
 *      Author: scheres
 */

#include "src/backprojector.h"

#ifdef TIMING
	#define RCTIC(timer,label) (timer.tic(label))
    #define RCTOC(timer,label) (timer.toc(label))
#else
	#define RCTIC(timer,label)
    #define RCTOC(timer,label)
#endif



void BackProjector::initialiseDataAndWeight(int current_size)
{

	initialiseData(current_size);
	weight.resize(data);

}

void BackProjector::initZeros(int current_size)
{

	initialiseDataAndWeight(current_size);
	data.initZeros();
	weight.initZeros();
}

void BackProjector::backproject2Dto3D(const MultidimArray<Complex > &f2d,
  	                              const Matrix2D<RFLOAT> &A,
                                      const MultidimArray<RFLOAT> *Mweight,
                                      RFLOAT r_ewald_sphere, bool is_positive_curvature,
                                      Matrix2D<RFLOAT>* magMatrix)
{
	RFLOAT m00, m10, m01, m11;

	if (magMatrix != 0)
	{
		m00 = (*magMatrix)(0,0);
		m10 = (*magMatrix)(1,0);
		m01 = (*magMatrix)(0,1);
		m11 = (*magMatrix)(1,1);
	}
	else
	{
		m00 = 1.0;
		m10 = 0.0;
		m01 = 0.0;
		m11 = 1.0;
	}

	// Use the inverse matrix
	Matrix2D<RFLOAT> Ainv;
	Ainv = A.inv();

	// Go from the 2D slice coordinates to the 3D coordinates
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	// max_r2 and min_r2_nn are defined in 3D-space
	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	const int min_r2_nn = ROUND(r_min_nn * padding_factor) * ROUND(r_min_nn * padding_factor);

	// precalculated coefficients for ellipse determination (see further down)

	// first, make sure A contains 2D distortion (lowercase 2D, uppercase 3D):
	const RFLOAT Am_Xx = Ainv(0,0) * m00 + Ainv(0,1) * m10;
	const RFLOAT Am_Xy = Ainv(0,0) * m01 + Ainv(0,1) * m11;
	const RFLOAT Am_Yx = Ainv(1,0) * m00 + Ainv(1,1) * m10;
	const RFLOAT Am_Yy = Ainv(1,0) * m01 + Ainv(1,1) * m11;
	const RFLOAT Am_Zx = Ainv(2,0) * m00 + Ainv(2,1) * m10;
	const RFLOAT Am_Zy = Ainv(2,0) * m01 + Ainv(2,1) * m11;

	// next, precompute (Am)^t Am into AtA:
	const RFLOAT AtA_xx = Am_Xx * Am_Xx + Am_Yx * Am_Yx + Am_Zx * Am_Zx;
	const RFLOAT AtA_xy = Am_Xx * Am_Xy + Am_Yx * Am_Yy + Am_Zx * Am_Zy;
	const RFLOAT AtA_yy = Am_Xy * Am_Xy + Am_Yy * Am_Yy + Am_Zy * Am_Zy;
	const RFLOAT AtA_xy2 = AtA_xy * AtA_xy;

//#define DEBUG_BACKP
#ifdef DEBUG_BACKP
	std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
	std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
	std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
	std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
	std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
	std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
	std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
	std::cerr << " r_max= "<< r_max << std::endl;
	std::cerr << " Ainv= " << Ainv << std::endl;
#endif

	// precalculate inverse of Ewald sphere diameter
	RFLOAT inv_diam_ewald = (r_ewald_sphere > 0.0)? 1.0 / (2.0 * r_ewald_sphere) : 0.0;

	if (!is_positive_curvature)
	{
		inv_diam_ewald *= -1.0;
	}

	const int s  = YSIZE(f2d);
	const int sh = XSIZE(f2d);

	for (int i = 0; i < s; i++)
	{
		int y, first_allowed_x;

		if (i < sh)
		{
			y = i;
			first_allowed_x = 0;
		}
		else
		{
			y = i - s;
			// x == 0 plane is stored twice in the FFTW format. Don't set it twice in backprojection!
			first_allowed_x = 1;
		}

		// Only iterate over the ellipse in the 2D-image corresponding to the sphere in 3D.
		// Find the x-range inside that ellipse for every given y:
		// |A*v|^2 <= R^2    (for v = (x,y)^t)
		// = v^t A^t A v =: v^t AtA v
		//   <=>
		// (AtA_xx) x^2 + (2 AtA_xy y) x + (AtA_yy y^2 - R^2) <= 0   (quadratic eq. in x)
		//   <=>
		// x in [q - d, q + d],
		// where: q := -AtA_xy y / AtA_xx,
		//        d := sqrt((AtA_xy y)^2 - AtA_xx (AtA_yy y^2 - R^2)) / AtA_xx

		RFLOAT discr = AtA_xy2 * y * y - AtA_xx * (AtA_yy * y * y - max_r2);

		if (discr < 0.0) continue; // no points inside ellipse for this y

		RFLOAT d = sqrt(discr) / AtA_xx;
		RFLOAT q = - AtA_xy * y / AtA_xx;

		int first_x = CEIL(q - d);
		int last_x = FLOOR(q + d);

		if (first_x < first_allowed_x) first_x = first_allowed_x;
		if (last_x > sh - 1) last_x = sh - 1;

		for (int x = first_x; x <= last_x; x++)
		{
			// Get the value from the input image
			Complex my_val = DIRECT_A2D_ELEM(f2d, i, x);

			RFLOAT my_weight;

			// Get the weight
			if (Mweight != NULL)
			{
				my_weight = DIRECT_A2D_ELEM(*Mweight, i, x);
			}
			else
			{
				my_weight = 1.0;
			}

			if (my_weight <= 0.) continue;

			/*
			In our implementation, (x, y) are not scaled because:

			x_on_ewald = x * r / sqrt(x * x + y * y + r * r)
			           = x / sqrt(1 + (x * x + y * y) / (r * r))
			           ~ x * (1 - (x * x + y * y) / (2 * r * r) + O(1/r^4)) # binomial expansion
			           = x + O(1/r^2)

			same for y_on_ewald

			z_on_ewald = r - r * r / sqrt(x * x + y * y + r * r)
			           ~ r - r * (1 - (x * x + y * y) / (2 * r * r) + O(1/r^4)) # binomial expansion
			           = (x * x + y * y) / (2 * r) + O(1/r^3)

			The error is < 0.0005 reciprocal voxel even for extreme cases
			like 200kV, 1500 A particle, 1 A / pix.
			*/

			// Get logical coordinates in the 3D map.
			// Make sure that the Ewald sphere is spherical even under anisotropic mag
			// by first undistorting (x,y) to obtain the true frequencies (xu,yu)

			RFLOAT xu = m00 * x + m01 * y;
			RFLOAT yu = m10 * x + m11 * y;

			RFLOAT z_on_ewaldp = inv_diam_ewald * (xu * xu + yu * yu);

			RFLOAT xp = Ainv(0,0) * xu + Ainv(0,1) * yu + Ainv(0,2) * z_on_ewaldp;
			RFLOAT yp = Ainv(1,0) * xu + Ainv(1,1) * yu + Ainv(1,2) * z_on_ewaldp;
			RFLOAT zp = Ainv(2,0) * xu + Ainv(2,1) * yu + Ainv(2,2) * z_on_ewaldp;

			double r2_3D = xp*xp + yp*yp + zp*zp;

			// redundant:
			if (r2_3D > max_r2)
			{
				continue;
			}

			if (interpolator == TRILINEAR || r2_3D < min_r2_nn)
			{
				bool is_neg_x;

				// Only asymmetric half is stored
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					zp = -zp;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
				int x0 = FLOOR(xp);
				RFLOAT fx = xp - x0;
				int x1 = x0 + 1;

				int y0 = FLOOR(yp);
				RFLOAT fy = yp - y0;
				y0 -=  STARTINGY(data);
				int y1 = y0 + 1;

				int z0 = FLOOR(zp);
				RFLOAT fz = zp - z0;
				z0 -= STARTINGZ(data);
				int z1 = z0 + 1;

				if (x0 < 0 || x0+1 >= data.xdim
				 || y0 < 0 || y0+1 >= data.ydim
				 || z0 < 0 || z0+1 >= data.zdim)
				{
					continue;
				}

				RFLOAT mfx = 1. - fx;
				RFLOAT mfy = 1. - fy;
				RFLOAT mfz = 1. - fz;

				RFLOAT dd000 = mfz * mfy * mfx;
				RFLOAT dd001 = mfz * mfy *  fx;
				RFLOAT dd010 = mfz *  fy * mfx;
				RFLOAT dd011 = mfz *  fy *  fx;
				RFLOAT dd100 =  fz * mfy * mfx;
				RFLOAT dd101 =  fz * mfy *  fx;
				RFLOAT dd110 =  fz *  fy * mfx;
				RFLOAT dd111 =  fz *  fy *  fx;

				if (is_neg_x)
				{
					my_val = conj(my_val);
				}

				// Store slice in 3D weighted sum
				DIRECT_A3D_ELEM(data, z0, y0, x0) += dd000 * my_val;
				DIRECT_A3D_ELEM(data, z0, y0, x1) += dd001 * my_val;
				DIRECT_A3D_ELEM(data, z0, y1, x0) += dd010 * my_val;
				DIRECT_A3D_ELEM(data, z0, y1, x1) += dd011 * my_val;
				DIRECT_A3D_ELEM(data, z1, y0, x0) += dd100 * my_val;
				DIRECT_A3D_ELEM(data, z1, y0, x1) += dd101 * my_val;
				DIRECT_A3D_ELEM(data, z1, y1, x0) += dd110 * my_val;
				DIRECT_A3D_ELEM(data, z1, y1, x1) += dd111 * my_val;
				// Store corresponding weights
				DIRECT_A3D_ELEM(weight, z0, y0, x0) += dd000 * my_weight;
				DIRECT_A3D_ELEM(weight, z0, y0, x1) += dd001 * my_weight;
				DIRECT_A3D_ELEM(weight, z0, y1, x0) += dd010 * my_weight;
				DIRECT_A3D_ELEM(weight, z0, y1, x1) += dd011 * my_weight;
				DIRECT_A3D_ELEM(weight, z1, y0, x0) += dd100 * my_weight;
				DIRECT_A3D_ELEM(weight, z1, y0, x1) += dd101 * my_weight;
				DIRECT_A3D_ELEM(weight, z1, y1, x0) += dd110 * my_weight;
				DIRECT_A3D_ELEM(weight, z1, y1, x1) += dd111 * my_weight;

			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR )
			{
				int x0 = ROUND(xp);
				int y0 = ROUND(yp);
				int z0 = ROUND(zp);

				bool is_neg_x;

				if (x0 < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					x0 = -x0;
					y0 = -y0;
					z0 = -z0;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

				const int xr = x0 - STARTINGX(data);
				const int yr = y0 - STARTINGY(data);
				const int zr = z0 - STARTINGZ(data);

				if (xr < 0 || xr >= data.xdim
				 || yr < 0 || yr >= data.ydim
				 || zr < 0 || zr >= data.zdim)
				{
					continue;
				}

				if (is_neg_x)
				{
					DIRECT_A3D_ELEM(data, zr, yr, xr) += conj(my_val);
					DIRECT_A3D_ELEM(weight, zr, yr, xr) += my_weight;
				}
				else
				{
					DIRECT_A3D_ELEM(data, zr, yr, xr) += my_val;
					DIRECT_A3D_ELEM(weight, zr, yr, xr) += my_weight;
				}
			} // endif NEAREST_NEIGHBOUR
			else
			{
				REPORT_ERROR("FourierInterpolator::backproject%%ERROR: unrecognized interpolator ");
			}
		} // endif x-loop
	} // endif y-loop
}

void BackProjector::backproject1Dto2D(const MultidimArray<Complex > &f1d,
                                      const Matrix2D<RFLOAT> &A,
                                      const MultidimArray<RFLOAT> *Mweight)
{
	Matrix2D<RFLOAT> Ainv = A.inv();
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	const int r_max_src = XSIZE(f1d) - 1;

	const int r_max_ref = r_max * padding_factor;
	const int r_max_ref_2 = r_max_ref * r_max_ref;

	// currently not used for some reason
	//const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

	for (int x = 0; x <= r_max_src; x++)
	{
		RFLOAT my_weight;

		if (Mweight != NULL)
		{
			my_weight = DIRECT_A1D_ELEM(*Mweight, x);

			if (my_weight <= 0.) continue;
		}
		else
		{
			my_weight = 1.;
		}

		Complex my_val = DIRECT_A1D_ELEM(f1d, x);

		// Get logical coordinates in the 3D map
		RFLOAT xp = Ainv(0,0) * x;
		RFLOAT yp = Ainv(1,0) * x;

		const RFLOAT r_ref_2 = xp*xp + yp*yp;

		if (r_ref_2 > r_max_ref_2) continue;

		if (interpolator == TRILINEAR /* && r_ref_2 < r_min_NN_ref_2*/)
		{
			// Only asymmetric half is stored
			const bool is_neg_x = xp < 0;

			if (is_neg_x)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
			}

			// Trilinear interpolation (with physical coords)
			// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
			// In that way use DIRECT_A2D_ELEM, rather than A2D_ELEM
			const int x0 = FLOOR(xp);
			const RFLOAT fx = xp - x0;
			const int x1 = x0 + 1;

			int y0 = FLOOR(yp);
			const RFLOAT fy = yp - y0;
			y0 -=  STARTINGY(data);
			const int y1 = y0 + 1;

			const RFLOAT mfx = 1. - fx;
			const RFLOAT mfy = 1. - fy;

			const RFLOAT dd00 = mfy * mfx;
			const RFLOAT dd01 = mfy *  fx;
			const RFLOAT dd10 =  fy * mfx;
			const RFLOAT dd11 =  fy *  fx;

			if (is_neg_x)
			{
				my_val = conj(my_val);
			}

			// Store slice in 3D weighted sum
			DIRECT_A2D_ELEM(data, y0, x0) += dd00 * my_val;
			DIRECT_A2D_ELEM(data, y0, x1) += dd01 * my_val;
			DIRECT_A2D_ELEM(data, y1, x0) += dd10 * my_val;
			DIRECT_A2D_ELEM(data, y1, x1) += dd11 * my_val;

			// Store corresponding weights
			DIRECT_A2D_ELEM(weight, y0, x0) += dd00 * my_weight;
			DIRECT_A2D_ELEM(weight, y0, x1) += dd01 * my_weight;
			DIRECT_A2D_ELEM(weight, y1, x0) += dd10 * my_weight;
			DIRECT_A2D_ELEM(weight, y1, x1) += dd11 * my_weight;

		} // endif TRILINEAR
		else if (interpolator == NEAREST_NEIGHBOUR )
		{
			const int x0 = ROUND(xp);
			const int y0 = ROUND(yp);

			if (x0 < 0)
			{
				A2D_ELEM(data, -y0, -x0) += conj(my_val);
				A2D_ELEM(weight, -y0, -x0) += my_weight;
			}
			else
			{
				A2D_ELEM(data, y0, x0) += my_val;
				A2D_ELEM(weight, y0, x0) += my_weight;
			}
		} // endif NEAREST_NEIGHBOUR
		else
		{
			REPORT_ERROR("FourierInterpolator::backproject1Dto2D%%ERROR: unrecognized interpolator ");
		}
	} // endif x-loop
}

void BackProjector::backrotate2D(const MultidimArray<Complex > &f2d,
                                 const Matrix2D<RFLOAT> &A,
                                 const MultidimArray<RFLOAT> *Mweight,
                                 Matrix2D<RFLOAT>* magMatrix)
{
	Matrix2D<RFLOAT> Ainv = A.inv();
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	RFLOAT m00, m10, m01, m11;

	if (magMatrix != 0)
	{
		m00 = (*magMatrix)(0,0);
		m10 = (*magMatrix)(1,0);
		m01 = (*magMatrix)(0,1);
		m11 = (*magMatrix)(1,1);
	}
	else
	{
		m00 = 1.0;
		m10 = 0.0;
		m01 = 0.0;
		m11 = 1.0;
	}

	const int r_max_ref = r_max * padding_factor;
	const int r_max_ref_2 = r_max_ref * r_max_ref;

	int min_r2_nn = r_min_nn * r_min_nn * padding_factor * padding_factor;

	// precalculated coefficients for ellipse determination (see further down)

	// first, make sure A contains 2D distortion (lowercase 2D, uppercase 3D):
	const RFLOAT Am_Xx = Ainv(0,0) * m00 + Ainv(0,1) * m10;
	const RFLOAT Am_Xy = Ainv(0,0) * m01 + Ainv(0,1) * m11;
	const RFLOAT Am_Yx = Ainv(1,0) * m00 + Ainv(1,1) * m10;
	const RFLOAT Am_Yy = Ainv(1,0) * m01 + Ainv(1,1) * m11;

	// next, precompute (Am)^t Am into AtA:
	const RFLOAT AtA_xx = Am_Xx * Am_Xx + Am_Yx * Am_Yx;
	const RFLOAT AtA_xy = Am_Xx * Am_Xy + Am_Yx * Am_Yy;
	const RFLOAT AtA_yy = Am_Xy * Am_Xy + Am_Yy * Am_Yy;
	const RFLOAT AtA_xy2 = AtA_xy * AtA_xy;

//#define DEBUG_BACKROTATE
#ifdef DEBUG_BACKROTATE
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
#endif
	const int s  = YSIZE(f2d);
	const int sh = XSIZE(f2d);

	for (int i = 0; i < s; i++)
	{
		int y, first_allowed_x;

		if (i < sh)
		{
			y = i;
			first_allowed_x = 0;
		}
		else
		{
			y = i - s;
			// x == 0 plane is stored twice in the FFTW format. Don't set it twice in backprojection!
			first_allowed_x = 1;
		}

		// Only iterate over the ellipse in the 2D-image corresponding to the sphere in 3D.
		// Find the x-range inside that ellipse for every given y:
		// |A*v|^2 <= R^2    (for v = (x,y)^t)
		// = v^t A^t A v =: v^t AtA v
		//   <=>
		// (AtA_xx) x^2 + (2 AtA_xy y) x + (AtA_yy y^2 - R^2) <= 0   (quadratic eq. in x)
		//   <=>
		// x in [q - d, q + d],
		// where: q := -AtA_xy y / AtA_xx,
		//        d := sqrt((AtA_xy y)^2 - AtA_xx (AtA_yy y^2 - R^2)) / AtA_xx

		RFLOAT discr = AtA_xy2 * y*y - AtA_xx * (AtA_yy * y*y - r_max_ref_2);

		if (discr < 0.0) continue; // no points inside ellipse for this y

		RFLOAT d = sqrt(discr) / AtA_xx;
		RFLOAT q = - AtA_xy * y / AtA_xx;

		int first_x = CEIL(q - d);
		int last_x = FLOOR(q + d);

		if (first_x < first_allowed_x) first_x = first_allowed_x;
		if (last_x > sh - 1) last_x = sh - 1;

		for (int x = first_x; x <= last_x; x++)
		{
			RFLOAT my_weight;

			if (Mweight != NULL)
			{
				my_weight = DIRECT_A2D_ELEM(*Mweight, i, x);

				if (my_weight <= 0.f) continue;
			}
			else
			{
				my_weight = 1.;
			}

			// Get the relevant value in the input image
			Complex my_val = DIRECT_A2D_ELEM(f2d, i, x);

			// Get logical coordinates in the 3D map
			RFLOAT xu = m00 * x + m01 * y;
			RFLOAT yu = m10 * x + m11 * y;

			RFLOAT xp = Ainv(0,0) * xu + Ainv(0,1) * yu;
			RFLOAT yp = Ainv(1,0) * xu + Ainv(1,1) * yu;

			RFLOAT r_ref_2 = xp * xp + yp * yp;

			if (interpolator == TRILINEAR || r_ref_2 < min_r2_nn)
			{
				const bool is_neg_x = xp < 0;
				// Only asymmetric half is stored
				if (is_neg_x)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A2D_ELEM, rather than A2D_ELEM
				const int x0 = FLOOR(xp);
				const RFLOAT fx = xp - x0;
				const int x1 = x0 + 1;

				int y0 = FLOOR(yp);
				const RFLOAT fy = yp - y0;
				y0 -=  STARTINGY(data);
				const int y1 = y0 + 1;

				const RFLOAT mfx = 1. - fx;
				const RFLOAT mfy = 1. - fy;

				const RFLOAT dd00 = mfy * mfx;
				const RFLOAT dd01 = mfy *  fx;
				const RFLOAT dd10 =  fy * mfx;
				const RFLOAT dd11 =  fy *  fx;

				if (is_neg_x)
				{
					my_val = conj(my_val);
				}

				// Store slice in 3D weighted sum
				DIRECT_A2D_ELEM(data, y0, x0) += dd00 * my_val;
				DIRECT_A2D_ELEM(data, y0, x1) += dd01 * my_val;
				DIRECT_A2D_ELEM(data, y1, x0) += dd10 * my_val;
				DIRECT_A2D_ELEM(data, y1, x1) += dd11 * my_val;

				// Store corresponding weights
				DIRECT_A2D_ELEM(weight, y0, x0) += dd00 * my_weight;
				DIRECT_A2D_ELEM(weight, y0, x1) += dd01 * my_weight;
				DIRECT_A2D_ELEM(weight, y1, x0) += dd10 * my_weight;
				DIRECT_A2D_ELEM(weight, y1, x1) += dd11 * my_weight;

			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR )
			{
				const int x0 = ROUND(xp);
				const int y0 = ROUND(yp);

				if (x0 < 0)
				{
					A2D_ELEM(data, -y0, -x0) += conj(my_val);
					A2D_ELEM(weight, -y0, -x0) += my_weight;
				}
				else
				{
					A2D_ELEM(data, y0, x0) += my_val;
					A2D_ELEM(weight, y0, x0) += my_weight;
				}
			} // endif NEAREST_NEIGHBOUR
			else
			{
				REPORT_ERROR("FourierInterpolator::backrotate2D%%ERROR: unrecognized interpolator ");
			}
		} // endif x-loop
	} // endif y-loop
}

void BackProjector::backrotate3D(const MultidimArray<Complex > &f3d,
                                 const Matrix2D<RFLOAT> &A,
                                 const MultidimArray<RFLOAT> *Mweight)
{
	// f3d should already be in the right size (ori_size,orihalfdim)
	// AND the points outside max_r should already be zero.

	Matrix2D<RFLOAT> Ainv = A.inv();
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	const int r_max_src = XSIZE(f3d) - 1;
	const int r_max_src_2 = r_max_src * r_max_src;

	const int r_max_ref = r_max * padding_factor;
	const int r_max_ref_2 = r_max_ref * r_max_ref;

	const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

//#define DEBUG_BACKROTATE
#ifdef DEBUG_BACKROTATE
	std::cerr << " XSIZE(f3d)= "<< XSIZE(f3d) << std::endl;
	std::cerr << " YSIZE(f3d)= "<< YSIZE(f3d) << std::endl;
	std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
	std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
	std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
	std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
	std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
	std::cerr << " max_r= "<< r_max << std::endl;
	std::cerr << " Ainv= " << Ainv << std::endl;
#endif

	for (int k = 0; k < ZSIZE(f3d); k++)
	{
		int z, x_min;

		// Don't search beyond square with side max_r
		if (k <= r_max_src)
		{
			z = k;
			x_min = 0;
		}
		else
		{
			z = k - ZSIZE(f3d);
			/// TODO: still check this better in the 3D case!!!
			// x==0 (y,z)-plane is stored twice in the FFTW format. Don't set it twice in BACKPROJECTION!
			x_min = 1;
		}

		int z2 = z * z;

		for (int i = 0; i < YSIZE(f3d); i++)
		{
			int y = (i <= r_max_src)? i : i - YSIZE(f3d);
			int y2 = y * y;

			const RFLOAT yz2 = y2 + z2;

			// avoid negative square root
			if (yz2 > r_max_src_2) continue;

			const int x_max = FLOOR(sqrt(r_max_src_2 - yz2));

			for (int x = x_min; x <= x_max; x++)
			{
				// Get logical coordinates in the 3D map
				RFLOAT xp = Ainv(0,0) * x + Ainv(0,1) * y + Ainv(0,2) * z;
				RFLOAT yp = Ainv(1,0) * x + Ainv(1,1) * y + Ainv(1,2) * z;
				RFLOAT zp = Ainv(2,0) * x + Ainv(2,1) * y + Ainv(2,2) * z;

				const int r_ref_2 = xp*xp + yp*yp + zp*zp;

				if (r_ref_2 > r_max_ref_2) continue;

				RFLOAT my_weight;

				// Get the weight
				if (Mweight != NULL)
				{
					my_weight = DIRECT_A3D_ELEM(*Mweight, k, i, x);

					if (my_weight <= 0.) continue;
				}
				else
				{
					my_weight = 1.;
				}

				Complex my_val = DIRECT_A3D_ELEM(f3d, k, i, x);

				if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2)
				{
					// Only asymmetric half is stored
					bool is_neg_x = xp < 0;

					if (is_neg_x)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
					}

					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					const int x0 = FLOOR(xp);
					const RFLOAT fx = xp - x0;
					const int x1 = x0 + 1;

					int y0 = FLOOR(yp);
					const RFLOAT fy = yp - y0;
					y0 -=  STARTINGY(data);
					const int y1 = y0 + 1;

					int z0 = FLOOR(zp);
					const RFLOAT fz = zp - z0;
					z0 -=  STARTINGZ(data);
					const int z1 = z0 + 1;

					const RFLOAT mfx = 1. - fx;
					const RFLOAT mfy = 1. - fy;
					const RFLOAT mfz = 1. - fz;

					const RFLOAT dd000 = mfz * mfy * mfx;
					const RFLOAT dd001 = mfz * mfy *  fx;
					const RFLOAT dd010 = mfz *  fy * mfx;
					const RFLOAT dd011 = mfz *  fy *  fx;
					const RFLOAT dd100 =  fz * mfy * mfx;
					const RFLOAT dd101 =  fz * mfy *  fx;
					const RFLOAT dd110 =  fz *  fy * mfx;
					const RFLOAT dd111 =  fz *  fy *  fx;

					if (is_neg_x)
					{
						my_val = conj(my_val);
					}

					// Store slice in 3D weighted sum
					DIRECT_A3D_ELEM(data, z0, y0, x0) += dd000 * my_val;
					DIRECT_A3D_ELEM(data, z0, y0, x1) += dd001 * my_val;
					DIRECT_A3D_ELEM(data, z0, y1, x0) += dd010 * my_val;
					DIRECT_A3D_ELEM(data, z0, y1, x1) += dd011 * my_val;
					DIRECT_A3D_ELEM(data, z1, y0, x0) += dd100 * my_val;
					DIRECT_A3D_ELEM(data, z1, y0, x1) += dd101 * my_val;
					DIRECT_A3D_ELEM(data, z1, y1, x0) += dd110 * my_val;
					DIRECT_A3D_ELEM(data, z1, y1, x1) += dd111 * my_val;

					// Store corresponding weights
					DIRECT_A3D_ELEM(weight, z0, y0, x0) += dd000 * my_weight;
					DIRECT_A3D_ELEM(weight, z0, y0, x1) += dd001 * my_weight;
					DIRECT_A3D_ELEM(weight, z0, y1, x0) += dd010 * my_weight;
					DIRECT_A3D_ELEM(weight, z0, y1, x1) += dd011 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y0, x0) += dd100 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y0, x1) += dd101 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y1, x0) += dd110 * my_weight;
					DIRECT_A3D_ELEM(weight, z1, y1, x1) += dd111 * my_weight;

				} // endif TRILINEAR
				else if (interpolator == NEAREST_NEIGHBOUR )
				{
					const int x0 = ROUND(xp);
					const int y0 = ROUND(yp);
					const int z0 = ROUND(zp);

					if (x0 < 0)
					{
						A3D_ELEM(data, -z0, -y0, -x0) += conj(my_val);
						A3D_ELEM(weight, -z0, -y0, -x0) += my_weight;
					}
					else
					{
						A3D_ELEM(data, z0, y0, x0) += my_val;
						A3D_ELEM(weight, z0, y0, x0) += my_weight;
					}

				} // endif NEAREST_NEIGHBOUR
				else
				{
					REPORT_ERROR("BackProjector::backrotate3D%%ERROR: unrecognized interpolator ");
				}
			} // endif x-loop
		} // endif y-loop
	} // endif z-loop
}

void BackProjector::getLowResDataAndWeight(MultidimArray<Complex > &lowres_data, MultidimArray<RFLOAT> &lowres_weight,
                                           int lowres_r_max)
{

	const int lowres_r2_max = ROUND(padding_factor * lowres_r_max) * ROUND(padding_factor * lowres_r_max);
	const int lowres_pad_size = 2 * (ROUND(padding_factor * lowres_r_max) + 1) + 1;

	// Check lowres_r_max is not too big
	if (lowres_r_max > r_max)
		REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");

	// Initialize lowres_data and low_res_weight arrays
	lowres_data.clear();
	lowres_weight.clear();

	if (ref_dim == 2)
	{
		lowres_data.resize(lowres_pad_size, lowres_pad_size / 2 + 1);
		lowres_weight.resize(lowres_pad_size, lowres_pad_size / 2 + 1);
	}
	else
	{
		lowres_data.resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
		lowres_weight.resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
	}
	lowres_data.setXmippOrigin();
	lowres_data.xinit=0;
	lowres_weight.setXmippOrigin();
	lowres_weight.xinit=0;

	// fill lowres arrays with relevant values
	FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data)
	{
		if (k*k + i*i + j*j <= lowres_r2_max)
		{
			A3D_ELEM(lowres_data, k, i, j) = A3D_ELEM(data, k , i, j);
			A3D_ELEM(lowres_weight, k, i, j) = A3D_ELEM(weight, k , i, j);
		}
	}
}

void BackProjector::setLowResDataAndWeight(MultidimArray<Complex > &lowres_data, MultidimArray<RFLOAT> &lowres_weight,
                                           int lowres_r_max)
{

	const int lowres_r2_max = ROUND(padding_factor * lowres_r_max) * ROUND(padding_factor * lowres_r_max);
	const int lowres_pad_size = 2 * (ROUND(padding_factor * lowres_r_max) + 1) + 1;

	// Check lowres_r_max is not too big
	if (lowres_r_max > r_max)
		REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");

	// Check sizes of lowres_data and lowres_weight
	if (YSIZE(lowres_data) != lowres_pad_size || XSIZE(lowres_data) != lowres_pad_size / 2 + 1 ||
			(ref_dim ==3 && ZSIZE(lowres_data) != lowres_pad_size) )
		REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_data is not of expected size...");
	if (YSIZE(lowres_weight) != lowres_pad_size || XSIZE(lowres_weight) != lowres_pad_size / 2 + 1 ||
			(ref_dim ==3 && ZSIZE(lowres_weight) != lowres_pad_size) )
		REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_weight is not of expected size...");

	// Re-set origin to the expected place
	lowres_data.setXmippOrigin();
	lowres_data.xinit=0;
	lowres_weight.setXmippOrigin();
	lowres_weight.xinit=0;

	// Overwrite data and weight with the lowres arrays
	FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data)
	{
		if (k*k + i*i + j*j <= lowres_r2_max)
		{
			A3D_ELEM(data, k, i, j) = A3D_ELEM(lowres_data, k , i, j);
			A3D_ELEM(weight, k, i, j) = A3D_ELEM(lowres_weight, k , i, j);
		}
	}
}

void BackProjector::getDownsampledAverage(MultidimArray<Complex>& avg, bool divide) const
{
	MultidimArray<RFLOAT> down_weight;

	// Pre-set down_data and down_weight sizes
	const int down_size = 2 * (r_max + 1) + 1;

	// Short side of data array
	switch (ref_dim)
	{
	case 2:
		avg.initZeros(down_size, down_size / 2 + 1);
		break;
	case 3:
		avg.initZeros(down_size, down_size, down_size / 2 + 1);
		break;
	default:
		REPORT_ERROR("BackProjector::getDownsampledAverage%%ERROR: Dimension of the data array should be 2 or 3");
	}
	// Set origin in the y.z-center, but on the left side for x.
	avg.setXmippOrigin();
	avg.xinit=0;
	// Resize down_weight the same as down_data
	down_weight.initZeros(avg);

	// Now calculate the down-sized sum
	int kp, ip, jp;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
	{
		kp = ROUND((RFLOAT)k/padding_factor);
		ip = ROUND((RFLOAT)i/padding_factor);
		jp = ROUND((RFLOAT)j/padding_factor);

// TMP
//#define CHECK_SIZE
#ifdef CHECK_SIZE
		if (kp > FINISHINGZ(avg) || ip > FINISHINGY(avg) || jp > FINISHINGX(avg) ||
		    kp < STARTINGZ(avg) || ip < STARTINGY(avg) || jp < STARTINGX(avg))
		{
			std::cerr << " kp= " << kp << " ip= " << ip << " jp= " << jp << std::endl;
			avg.printShape();
			REPORT_ERROR("BackProjector::getDownsampledAverage: indices out of range");
		}
#endif
		A3D_ELEM(avg, kp, ip, jp) += A3D_ELEM(data, k , i, j);
		A3D_ELEM(down_weight, kp, ip, jp) += (divide? A3D_ELEM(weight, k , i, j) : 1.0);
	}

	// Calculate the straightforward average in the downsampled arrays
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(avg)
	{
		if (DIRECT_MULTIDIM_ELEM(down_weight, n) > 0.)
		{
			DIRECT_MULTIDIM_ELEM(avg, n) /= DIRECT_MULTIDIM_ELEM(down_weight, n);
		}
		else
		{
			DIRECT_MULTIDIM_ELEM(avg, n) = 0.;
		}
	}
}

void BackProjector::calculateDownSampledFourierShellCorrelation(const MultidimArray<Complex>& avg1,
                                                                const MultidimArray<Complex>& avg2,
                                                                MultidimArray<RFLOAT>& fsc) const
{
	if (!avg1.sameShape(avg2))
		REPORT_ERROR("ERROR BackProjector::calculateDownSampledFourierShellCorrelation: two arrays have different sizes");

	MultidimArray<RFLOAT> num, den1, den2;
	num.initZeros(ori_size/2 + 1);
	den1.initZeros(num);
	den2.initZeros(num);
	fsc.initZeros(num);

	FOR_ALL_ELEMENTS_IN_ARRAY3D(avg1)
	{
		const RFLOAT R = sqrt(k*k + i*i + j*j);

		if (R > r_max) continue;

		int idx = ROUND(R);

		Complex z1 = A3D_ELEM(avg1, k, i, j);
		Complex z2 = A3D_ELEM(avg2, k, i, j);

		RFLOAT nrmz1 = z1.norm();
		RFLOAT nrmz2 = z2.norm();

		num(idx) += z1.real * z2.real + z1.imag * z2.imag;
		den1(idx) += nrmz1;
		den2(idx) += nrmz2;
	}

	FOR_ALL_ELEMENTS_IN_ARRAY1D(fsc)
	{
		if (den1(i)*den2(i) > 0.)
		{
			fsc(i) = num(i)/sqrt(den1(i)*den2(i));
		}
	}

	// Always set zero-resolution shell to FSC=1
	// Raimond Ravelli reported a problem with FSC=1 at res=0 on 13feb2013...
	// (because of a suboptimal normalisation scheme, but anyway)
	fsc(0) = 1.;
}

void BackProjector::updateSSNRarrays(RFLOAT tau2_fudge,
                                     MultidimArray<RFLOAT> &tau2_io,
                                     MultidimArray<RFLOAT> &sigma2_out,
                                     MultidimArray<RFLOAT> &data_vs_prior_out,
                                     MultidimArray<RFLOAT> &fourier_coverage_out,
                                     const MultidimArray<RFLOAT>& fsc,
									 const MultidimArray<RFLOAT>& avgctf2,
                                     bool update_tau2_with_fsc,
                                     bool is_whole_instead_of_half,
									 bool correct_tau2_by_avgctf2)
{
	// never rely on references (handed to you from the outside) for computation:
	// they could be the same (i.e. reconstruct(..., dummy, dummy, dummy, dummy, ...); )
	MultidimArray<RFLOAT> sigma2, data_vs_prior, fourier_coverage;
	MultidimArray<RFLOAT> tau2 = tau2_io;
	MultidimArray<RFLOAT> counter;
	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	RFLOAT oversampling_correction = (ref_dim == 3) ? (padding_factor * padding_factor * padding_factor) : (padding_factor * padding_factor);

	// First calculate the radial average of the (inverse of the) power of the noise in the reconstruction
	// This is the left-hand side term in the nominator of the Wiener-filter-like update formula
	// and it is stored inside the weight vector
	// Then, if (do_map) add the inverse of tau2-spectrum values to the weight
	sigma2.initZeros(ori_size/2 + 1);
	counter.initZeros(ori_size/2 + 1);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(weight)
	{
		const int r2 = k * k + i * i + j * j;
		if (r2 < max_r2)
		{
			int ires = ROUND(sqrt((RFLOAT)r2) / padding_factor);
			RFLOAT invw = oversampling_correction * A3D_ELEM(weight, k, i, j);
			DIRECT_A1D_ELEM(sigma2, ires) += invw;
			DIRECT_A1D_ELEM(counter, ires) += 1.;
		}
	}

	// Average (inverse of) sigma2 in reconstruction
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sigma2)
	{
		if (DIRECT_A1D_ELEM(sigma2, i) > 1e-10)
			DIRECT_A1D_ELEM(sigma2, i) = DIRECT_A1D_ELEM(counter, i) / DIRECT_A1D_ELEM(sigma2, i);
		else if (DIRECT_A1D_ELEM(sigma2, i) == 0)
			DIRECT_A1D_ELEM(sigma2, i) = 0.;
		else
		{
			std::cerr << " DIRECT_A1D_ELEM(sigma2, i)= " << DIRECT_A1D_ELEM(sigma2, i) << std::endl;
			REPORT_ERROR("BackProjector::reconstruct: ERROR: unexpectedly small, yet non-zero sigma2 value, this should not happen...");
		}
	}

	tau2.reshape(ori_size/2 + 1);
	data_vs_prior.initZeros(ori_size/2 + 1);
	fourier_coverage.initZeros(ori_size/2 + 1);
	counter.initZeros(ori_size/2 + 1);
	if (update_tau2_with_fsc)
	{
		// Then calculate new tau2 values, based on the FSC
		if (!fsc.sameShape(sigma2) || !fsc.sameShape(tau2))
		{
			fsc.printShape(std::cerr);
			tau2.printShape(std::cerr);
			sigma2.printShape(std::cerr);
			REPORT_ERROR("ERROR BackProjector::reconstruct: sigma2, tau2 and fsc have different sizes");
		}
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sigma2)
		{
			// FSC cannot be negative or zero for conversion into tau2
			RFLOAT myfsc = XMIPP_MAX(0.001, DIRECT_A1D_ELEM(fsc, i));
			if (is_whole_instead_of_half)
			{
				// Factor two because of twice as many particles
				// Sqrt-term to get 60-degree phase errors....
				myfsc = sqrt(2. * myfsc / (myfsc + 1.));
			}
			myfsc = XMIPP_MIN(0.999, myfsc);
			RFLOAT myssnr = myfsc / (1. - myfsc);
			// Sjors 29nov2017 try tau2_fudge for pulling harder on Refine3D runs...
			myssnr *= tau2_fudge;
			RFLOAT fsc_based_tau = myssnr * DIRECT_A1D_ELEM(sigma2, i);
			DIRECT_A1D_ELEM(tau2, i) = fsc_based_tau;
			// data_vs_prior is merely for reporting: it is not used for anything in the reconstruction
			DIRECT_A1D_ELEM(data_vs_prior, i) = myssnr;
		}
	}

	// Now accumulate data_vs_prior if (!update_tau2_with_fsc)
	// Also accumulate fourier_coverage
	FOR_ALL_ELEMENTS_IN_ARRAY3D(weight)
	{
		int r2 = k * k + i * i + j * j;
		if (r2 < max_r2)
		{
			int ires = ROUND( sqrt((RFLOAT)r2) / padding_factor );
			RFLOAT invw = A3D_ELEM(weight, k, i, j);

			RFLOAT invtau2;
			if (DIRECT_A1D_ELEM(tau2, ires) > 0.)
			{
				// Calculate inverse of tau2
				invtau2 = 1. / (oversampling_correction * tau2_fudge * DIRECT_A1D_ELEM(tau2, ires));

				// Correct tau2 estimate by avgctf2 (for ctf_premultiplied images)
				if (correct_tau2_by_avgctf2 && DIRECT_A1D_ELEM(avgctf2, ires) > 0.)
					invtau2 *= 1. / DIRECT_A1D_ELEM(avgctf2, ires);

			}
			else if (DIRECT_A1D_ELEM(tau2, ires) == 0.)
			{
				// If tau2 is zero, use small value instead
				invtau2 = 1./ ( 0.001 * invw);
			}
			else
			{
				std::cerr << " sigma2= " << sigma2 << std::endl;
				std::cerr << " fsc= " << fsc << std::endl;
				std::cerr << " tau2= " << tau2 << std::endl;
				REPORT_ERROR("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
			}

			// Keep track of spectral evidence-to-prior ratio and remaining noise in the reconstruction
			if (!update_tau2_with_fsc)
			{
				DIRECT_A1D_ELEM(data_vs_prior, ires) += invw / invtau2;
			}

			// Keep track of the coverage in Fourier space
			if (invw / invtau2 >= 1.)
			{
				DIRECT_A1D_ELEM(fourier_coverage, ires) += 1.;
			}

			DIRECT_A1D_ELEM(counter, ires) += 1.;

		}
	}

	// Average data_vs_prior
	if (!update_tau2_with_fsc)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(data_vs_prior)
		{
			if (i > r_max)
				DIRECT_A1D_ELEM(data_vs_prior, i) = 0.;
			else if (DIRECT_A1D_ELEM(counter, i) < 0.001)
				DIRECT_A1D_ELEM(data_vs_prior, i) = 999.;
			else
				DIRECT_A1D_ELEM(data_vs_prior, i) /= DIRECT_A1D_ELEM(counter, i);
		}
	}

	// Calculate Fourier coverage in each shell
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fourier_coverage)
	{
		if (DIRECT_A1D_ELEM(counter, i) > 0.)
			DIRECT_A1D_ELEM(fourier_coverage, i) /= DIRECT_A1D_ELEM(counter, i);
	}

	// Send back the output
	tau2_io = tau2;
	sigma2_out = sigma2;
	data_vs_prior_out = data_vs_prior;
	fourier_coverage_out = fourier_coverage;
}

void BackProjector::externalReconstruct(MultidimArray<RFLOAT> &vol_out,
                                        FileName &fn_out,
                                        MultidimArray<RFLOAT> &fsc_halves_io,
                                        MultidimArray<RFLOAT> &tau2_io,
										MultidimArray<RFLOAT> &sigma2_ref,
										MultidimArray<RFLOAT> &data_vs_prior,
										bool is_whole_instead_of_half,
                                        RFLOAT tau2_fudge,
                                        int verb)
{

	FileName fn_recons = fn_out+"_external_reconstruct.mrc";
	FileName fn_star = fn_out+"_external_reconstruct.star";
	FileName fn_out_star = fn_out+"_external_reconstruct_out.star";
	MultidimArray<RFLOAT> fsc_halves = fsc_halves_io;
	MultidimArray<RFLOAT> tau2 = tau2_io;

	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	int padoridim = ROUND(padding_factor * ori_size);

	// Write out data array
	Image<Complex> Idata;
	if (ref_dim == 2) Idata().resize(pad_size, pad_size/2+1);
	else Idata().resize(pad_size, pad_size, pad_size/2+1);
	Projector::decenter(data, Idata(), max_r2);
	windowFourierTransform(Idata(), padoridim);
	ComplexIO::write(Idata(), fn_out+"_external_reconstruct_data", ".mrc");
	Idata.clear();

	// Write out weight array
	Image<RFLOAT> Iweight;
	if (ref_dim == 2) Iweight().resize(pad_size, pad_size/2+1);
	else Iweight().resize(pad_size, pad_size, pad_size/2+1);
	Projector::decenter(weight, Iweight(), max_r2);
	windowFourierTransform(Iweight(), padoridim);
	Iweight.write(fn_out+"_external_reconstruct_weight.mrc");
	Iweight.clear();

	// Write out STAR file for input to external reconstruction program
	MetaDataTable MDlist, MDtau;

	MDlist.setName("external_reconstruct_general");
	MDlist.setIsList(true);
	MDlist.addObject();
	MDlist.setValue(EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_REAL, fn_out+"_external_reconstruct_data_real.mrc");
	MDlist.setValue(EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_IMAG, fn_out+"_external_reconstruct_data_imag.mrc");
	MDlist.setValue(EMDL_OPTIMISER_EXTERNAL_RECONS_WEIGHT, fn_out+"_external_reconstruct_weight.mrc");
	MDlist.setValue(EMDL_OPTIMISER_EXTERNAL_RECONS_RESULT, fn_recons);
	MDlist.setValue(EMDL_OPTIMISER_EXTERNAL_RECONS_NEWSTAR, fn_out_star);
	MDlist.setValue(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge);
	MDlist.setValue(EMDL_MLMODEL_PADDING_FACTOR, padding_factor);
	MDlist.setValue(EMDL_MLMODEL_DIMENSIONALITY, ref_dim);
	MDlist.setValue(EMDL_MLMODEL_ORIGINAL_SIZE, ori_size);
	MDlist.setValue(EMDL_MLMODEL_CURRENT_SIZE, 2*r_max);

	MDtau.setName("external_reconstruct_tau2");
	for (int ii = 0; ii < XSIZE(tau2); ii++)
	{
		MDtau.addObject();
		MDtau.setValue(EMDL_SPECTRAL_IDX, ii);
		MDtau.setValue(EMDL_MLMODEL_TAU2_REF, tau2(ii));
		MDtau.setValue(EMDL_MLMODEL_FSC_HALVES_REF, fsc_halves(ii));
	}

	std::ofstream  fh;
	fh.open((fn_star).c_str(), std::ios::out);
	if (!fh) REPORT_ERROR( (std::string)"BackProjector::externalReconstruct: Cannot write file: " + fn_star);
	MDlist.write(fh);
	MDtau.write(fh);
	fh.close();


	// Make the system call: program name plus the STAR file for the external reconstruction program as its first argument
	char *my_exec = getenv ("RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE");
	char default_exec[]=DEFAULT_EXTERNAL_RECONSTRUCT;
	if (my_exec == NULL)
	{
		my_exec = default_exec;
	}
	std::string command = std::string(my_exec) + " " + fn_star;

	if (verb > 0) std::cout << std::endl << " + Making system call for external reconstruction: " << command << std::endl;

	int res = system(command.c_str());
	if (res) REPORT_ERROR(" ERROR: there was something wrong with system call: " + command);
	else if (verb > 0) std::cout << " + External reconstruction finished successfully, reading result back in ... " << std::endl;

	// Read the resulting map back into memory
	Iweight.read(fn_recons);
	vol_out = Iweight();
	vol_out.setXmippOrigin();

	if (exists(fn_out_star))
	{
		MetaDataTable MDnewtau;
		MDnewtau.read(fn_out_star);
		if (!MDnewtau.containsLabel(EMDL_SPECTRAL_IDX))
			REPORT_ERROR("ERROR: external reconstruct output STAR file does not contain spectral idx!");

		// Directly update tau2 spectrum
		int idx;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDnewtau)
		{
			MDnewtau.getValue(EMDL_SPECTRAL_IDX, idx);
			if (idx >= XSIZE(tau2_io)) continue;
			if (MDnewtau.containsLabel(EMDL_MLMODEL_TAU2_REF))
			{
				MDnewtau.getValue(EMDL_MLMODEL_TAU2_REF, tau2_io(idx));
				data_vs_prior(idx) = tau2_io(idx) / sigma2_ref(idx);
			}
			else if (MDnewtau.containsLabel(EMDL_POSTPROCESS_FSC_GENERAL))
			{
				MDnewtau.getValue(EMDL_SPECTRAL_IDX, idx);
				MDnewtau.getValue(EMDL_POSTPROCESS_FSC_GENERAL, fsc_halves_io(idx));

				RFLOAT myfsc = XMIPP_MAX(0.001, fsc_halves_io(idx));
				if (is_whole_instead_of_half)
				{
					// Factor two because of twice as many particles
					// Sqrt-term to get 60-degree phase errors....
					myfsc = sqrt(2. * myfsc / (myfsc + 1.));
				}
				myfsc = XMIPP_MIN(0.999, myfsc);
				RFLOAT myssnr = myfsc / (1. - myfsc);
				myssnr *= tau2_fudge;
				tau2_io(idx) = myssnr * sigma2_ref(idx);
				data_vs_prior(idx) = myssnr;
			}
			else
			{
				REPORT_ERROR("ERROR: output STAR file from external reconstruct does not contain tau2 or FSC array");
			}
		}
		if (verb > 0) std::cout << " + External reconstruction successfully updated external tau2 array ... " << std::endl;
	}

}

void BackProjector::reconstruct(MultidimArray<RFLOAT> &vol_out,
                                int max_iter_preweight,
                                bool do_map,
                                const MultidimArray<RFLOAT> &tau2,
                                RFLOAT tau2_fudge,
                                RFLOAT normalise,
                                int minres_map,
                                bool printTimes,
                                Image<RFLOAT>* weight_out)
{
#ifdef TIMING
	Timer ReconTimer;
	int ReconS_1 = ReconTimer.setNew(" RcS1_Init ");
	int ReconS_2 = ReconTimer.setNew(" RcS2_Shape&Noise ");
	int ReconS_2_5 = ReconTimer.setNew(" RcS2.5_Regularize ");
	int ReconS_3 = ReconTimer.setNew(" RcS3_skipGridding ");
	int ReconS_4 = ReconTimer.setNew(" RcS4_doGridding_norm ");
	int ReconS_5 = ReconTimer.setNew(" RcS5_doGridding_init ");
	int ReconS_6 = ReconTimer.setNew(" RcS6_doGridding_iter ");
	int ReconS_7 = ReconTimer.setNew(" RcS7_doGridding_apply ");
	int ReconS_8 = ReconTimer.setNew(" RcS8_blobConvolute ");
	int ReconS_9 = ReconTimer.setNew(" RcS9_blobResize ");
	int ReconS_10 = ReconTimer.setNew(" RcS10_blobSetReal ");
	int ReconS_11 = ReconTimer.setNew(" RcS11_blobSetTemp ");
	int ReconS_12 = ReconTimer.setNew(" RcS12_blobTransform ");
	int ReconS_13 = ReconTimer.setNew(" RcS13_blobCenterFFT ");
	int ReconS_14 = ReconTimer.setNew(" RcS14_blobNorm1 ");
	int ReconS_15 = ReconTimer.setNew(" RcS15_blobSoftMask ");
	int ReconS_16 = ReconTimer.setNew(" RcS16_blobNorm2 ");
	int ReconS_17 = ReconTimer.setNew(" RcS17_WindowReal ");
	int ReconS_18 = ReconTimer.setNew(" RcS18_GriddingCorrect ");
	int ReconS_19 = ReconTimer.setNew(" RcS19_tauInit ");
	int ReconS_20 = ReconTimer.setNew(" RcS20_tausetReal ");
	int ReconS_21 = ReconTimer.setNew(" RcS21_tauTransform ");
	int ReconS_22 = ReconTimer.setNew(" RcS22_tautauRest ");
	int ReconS_23 = ReconTimer.setNew(" RcS23_tauShrinkToFit ");
	int ReconS_24 = ReconTimer.setNew(" RcS24_extra ");
#endif

	RCTIC(ReconTimer,ReconS_1);
	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	RFLOAT oversampling_correction = (ref_dim == 3) ? (padding_factor * padding_factor * padding_factor) : (padding_factor * padding_factor);

//#define DEBUG_RECONSTRUCT
#ifdef DEBUG_RECONSTRUCT
	Image<RFLOAT> ttt;
	FileName fnttt;
	ttt()=weight;
	ttt.write("reconstruct_initial_weight.spi");
	std::cerr << " pad_size= " << pad_size << " padding_factor= " << padding_factor << " max_r2= " << max_r2 << std::endl;
#endif

	// Set Fconv to the right size
	if (ref_dim == 2)
		vol_out.setDimensions(pad_size, pad_size, 1, 1);
	else
        // Too costly to actually allocate the space
        // Trick transformer with the right dimensions
        vol_out.setDimensions(pad_size, pad_size, pad_size, 1);

	FourierTransformer transformer;
	transformer.setReal(vol_out); // Fake set real. 1. Allocate space for Fconv 2. calculate plans.
	MultidimArray<Complex>& Fconv = transformer.getFourierReference();
	vol_out.clear(); // Reset dimensions to 0

	RCTOC(ReconTimer,ReconS_1);
	RCTIC(ReconTimer,ReconS_2);

	// Go from projector-centered to FFTW-uncentered
	MultidimArray<RFLOAT> Fweight;
	Fweight.reshape(Fconv);
	Projector::decenter(weight, Fweight, max_r2);

	RCTOC(ReconTimer,ReconS_2);
	RCTIC(ReconTimer,ReconS_2_5);
	// Apply MAP-additional term to the Fweight array
	// This will regularise the actual reconstruction
	if (do_map)
	{
		// Then, add the inverse of tau2-spectrum values to the weight
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
 		{
			int r2 = kp * kp + ip * ip + jp * jp;
			if (r2 < max_r2)
			{
				int ires = ROUND(sqrt((RFLOAT)r2) / padding_factor);
				RFLOAT invw = DIRECT_A3D_ELEM(Fweight, k, i, j);

				RFLOAT invtau2;
				if (DIRECT_A1D_ELEM(tau2, ires) > 0.)
				{
					// Calculate inverse of tau2
					invtau2 = 1. / (oversampling_correction * tau2_fudge * DIRECT_A1D_ELEM(tau2, ires));
				}
				else if (DIRECT_A1D_ELEM(tau2, ires) < 1e-20)
				{
					// If tau2 is zero, use small value instead
					if (invw > 1e-20) invtau2 = 1./ ( 0.001 * invw);
					else invtau2 = 0.;
				}
				else
				{
					std::cerr << " tau2= " << tau2 << std::endl;
					REPORT_ERROR("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
				}

				// Only for (ires >= minres_map) add Wiener-filter like term
				if (ires >= minres_map)
				{
					// Now add the inverse-of-tau2_class term
					invw += invtau2;
					// Store the new weight again in Fweight
					DIRECT_A3D_ELEM(Fweight, k, i, j) = invw;
				}
			}
		}
	} //end if do_map

	RCTOC(ReconTimer,ReconS_2_5);
	if (skip_gridding || max_iter_preweight <= 0)
	{
		RCTIC(ReconTimer,ReconS_3);
		Fconv.initZeros(); // to remove any stuff from the input volume
		Projector::decenter(data, Fconv, max_r2);

		// Prevent divisions by zero: set Fweight to at least 1/1000th of the radially averaged weight at that resolution
		// beyond r_max, set Fweight to at least 1/1000th of the radially averaged weight at r_max;
		MultidimArray<RFLOAT> radavg_weight(r_max), counter(r_max);
		const int round_max_r2 = ROUND(r_max * padding_factor * r_max * padding_factor);
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fweight)
		{
			const int r2 = kp * kp + ip * ip + jp * jp;
			// Note that (r < ires) != (r2 < max_r2), because max_r2 = ROUND(r_max * padding_factor)^2.
			// We have to use round_max_r2 = ROUND((r_max * padding_factor)^2).
			// e.g. k = 0, i = 7, j = 28, max_r2 = 841, r_max = 16, padding_factor = 18.
			if (r2 < round_max_r2)
			{
				const int ires = FLOOR(sqrt((RFLOAT)r2) / padding_factor);
				if (ires >= XSIZE(radavg_weight))
				{
					std::cerr << " k= " << k << " i= " << i << " j= " << j << std::endl;
					std::cerr << " ires= " << ires << " XSIZE(radavg_weight)= " << XSIZE(radavg_weight) << std::endl;
					REPORT_ERROR("BUG: ires >=XSIZE(radavg_weight) ");
				}
				DIRECT_A1D_ELEM(radavg_weight, ires) += DIRECT_A3D_ELEM(Fweight, k, i, j);
				DIRECT_A1D_ELEM(counter, ires) += 1.;
			}
		}

		// Calculate 1/1000th of radial averaged weight
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(radavg_weight)
		{
			if (DIRECT_A1D_ELEM(counter, i) > 0. || DIRECT_A1D_ELEM(radavg_weight, i) > 0.)
			{
				DIRECT_A1D_ELEM(radavg_weight, i) /= 1000.* DIRECT_A1D_ELEM(counter, i);
			}
			else
			{
				std::cerr << " counter= " << counter << std::endl;
				std::cerr << " radavg_weight= " << radavg_weight << std::endl;
				REPORT_ERROR("BUG: zeros in counter or radavg_weight!");
			}
		}

		bool have_warned = false;
		// perform XMIPP_MAX on all weight elements, and do division of data/weight
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fweight)
		{
			const int r2 = kp * kp + ip * ip + jp * jp;
			const int ires = FLOOR(sqrt((RFLOAT)r2) / padding_factor);
			const RFLOAT weight =  XMIPP_MAX(DIRECT_A3D_ELEM(Fweight, k, i, j), DIRECT_A1D_ELEM(radavg_weight, (ires < r_max) ? ires : (r_max - 1)));
			if (weight == 0.)
			{
				if (!have_warned && abs(DIRECT_A3D_ELEM(Fconv, k, i, j)) > 0.)
				{
					std::cerr << " WARNING: ignoring divide by zero in skip_gridding: ires = " << ires << " kp = " << kp << " ip = " << ip << " jp = " << jp << std::endl;
                    std::cerr << " Fconv= " << DIRECT_A3D_ELEM(Fconv, k, i, j) << " Fweight= " << DIRECT_A3D_ELEM(Fweight, k, i, j) << " radavg_weight=" <<DIRECT_A1D_ELEM(radavg_weight, (ires < r_max) ? ires : (r_max - 1)) << std::endl;
					std::cerr << " max_r2 = " << max_r2 << " r_max = " << r_max << " padding_factor = " << padding_factor
					           << " ROUND(sqrt(max_r2)) = " << ROUND(sqrt(max_r2)) << " ROUND(r_max * padding_factor) = " << ROUND(r_max * padding_factor) << std::endl;
					have_warned = true;
				}
			}
			else
			{
				DIRECT_A3D_ELEM(Fconv, k, i, j) /= weight;
			}
		}
	}
	else
	{
		RCTIC(ReconTimer,ReconS_4);

		// Divide both data and Fweight by normalisation factor to prevent FFT's with very large values....
#ifdef DEBUG_RECONSTRUCT
		std::cerr << " normalise= " << normalise << std::endl;
#endif
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fweight)
		{
			DIRECT_MULTIDIM_ELEM(Fweight, n) /= normalise;
		}
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
		{
			DIRECT_MULTIDIM_ELEM(data, n) /= normalise;
		}

		RCTOC(ReconTimer,ReconS_4);
		RCTIC(ReconTimer,ReconS_5);

		// Initialise Fnewweight with 1's and 0's. (also see comments below)
		FOR_ALL_ELEMENTS_IN_ARRAY3D(weight)
		{
			if (k * k + i * i + j * j < max_r2)
				A3D_ELEM(weight, k, i, j) = 1.;
			else
				A3D_ELEM(weight, k, i, j) = 0.;
		}
		// Fnewweight can become too large for a float: always keep this one in double-precision
		MultidimArray<double> Fnewweight;
		Fnewweight.reshape(Fconv);
		decenter(weight, Fnewweight, max_r2);

		RCTOC(ReconTimer,ReconS_5);
		// Iterative algorithm as in  Eq. [14] in Pipe & Menon (1999)
		// or Eq. (4) in Matej (2001)
		for (int iter = 0; iter < max_iter_preweight; iter++)
		{
			//std::cout << "    iteration " << (iter+1) << "/" << max_iter_preweight << "\n";
			RCTIC(ReconTimer,ReconS_6);
			// Set Fnewweight * Fweight in the transformer
			// In Matej et al (2001), weights w_P^i are convoluted with the kernel,
			// and the initial w_P^0 are 1 at each sampling point
			// Here the initial weights are also 1 (see initialisation Fnewweight above),
			// but each "sampling point" counts "Fweight" times!
			// That is why Fnewweight is multiplied by Fweight prior to the convolution

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
			{
				DIRECT_MULTIDIM_ELEM(Fconv, n) = DIRECT_MULTIDIM_ELEM(Fnewweight, n) * DIRECT_MULTIDIM_ELEM(Fweight, n);
			}

			// convolute through Fourier-transform (as both grids are rectangular)
			// Note that convoluteRealSpace acts on the complex array inside the transformer
			convoluteBlobRealSpace(transformer, false);

			RFLOAT w, corr_min = LARGE_NUMBER, corr_max = -LARGE_NUMBER, corr_avg=0., corr_nn=0.;

			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
			{
				if (kp * kp + ip * ip + jp * jp < max_r2)
				{

					// Make sure no division by zero can occur....
					w = XMIPP_MAX(1e-6, abs(DIRECT_A3D_ELEM(Fconv, k, i, j)));
					// Monitor min, max and avg conv_weight
					corr_min = XMIPP_MIN(corr_min, w);
					corr_max = XMIPP_MAX(corr_max, w);
					corr_avg += w;
					corr_nn += 1.;
					// Apply division of Eq. [14] in Pipe & Menon (1999)
					DIRECT_A3D_ELEM(Fnewweight, k, i, j) /= w;
				}
			}

            RCTOC(ReconTimer,ReconS_6);

#ifdef DEBUG_RECONSTRUCT
			std::cerr << " PREWEIGHTING ITERATION: "<< iter + 1 << " OF " << max_iter_preweight << std::endl;
			// report of maximum and minimum values of current conv_weight
			std::cerr << " corr_avg= " << corr_avg / corr_nn << std::endl;
			std::cerr << " corr_min= " << corr_min << std::endl;
			std::cerr << " corr_max= " << corr_max << std::endl;
#endif
		}

		RCTIC(ReconTimer,ReconS_7);

#ifdef DEBUG_RECONSTRUCT
		Image<double> tttt;
		tttt()=Fnewweight;
		tttt.write("reconstruct_gridding_weight.spi");
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
		{
			DIRECT_MULTIDIM_ELEM(ttt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fconv, n));
		}
		ttt.write("reconstruct_gridding_correction_term.spi");
#endif


		// Clear memory
		Fweight.clear();

		// Note that Fnewweight now holds the approximation of the inverse of the weights on a regular grid

		// Now do the actual reconstruction with the data array
		// Apply the iteratively determined weight
		Fconv.initZeros(); // to remove any stuff from the input volume
		Projector::decenter(data, Fconv, max_r2);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
		{
#ifdef  RELION_SINGLE_PRECISION
			// Prevent numerical instabilities in single-precision reconstruction with very unevenly sampled orientations
			if (DIRECT_MULTIDIM_ELEM(Fnewweight, n) > 1e20)
				DIRECT_MULTIDIM_ELEM(Fnewweight, n) = 1e20;
#endif
			DIRECT_MULTIDIM_ELEM(Fconv, n) *= DIRECT_MULTIDIM_ELEM(Fnewweight, n);
		}

		// Clear memory
		Fnewweight.clear();

		RCTOC(ReconTimer,ReconS_7);
	} // end if !skip_gridding

// Gridding theory says one now has to interpolate the fine grid onto the coarse one using a blob kernel
// and then do the inverse transform and divide by the FT of the blob (i.e. do the gridding correction)
// In practice, this gives all types of artefacts (perhaps I never found the right implementation?!)
// Therefore, window the Fourier transform and then do the inverse transform
//#define RECONSTRUCT_CONVOLUTE_BLOB
#ifdef RECONSTRUCT_CONVOLUTE_BLOB

	// Apply the same blob-convolution as above to the data array
	// Mask real-space map beyond its original size to prevent aliasing in the downsampling step below
	RCTIC(ReconTimer,ReconS_8);
	convoluteBlobRealSpace(transformer, true);
	RCTOC(ReconTimer,ReconS_8);
	RCTIC(ReconTimer,ReconS_9);
	// Now just pick every 3rd pixel in Fourier-space (i.e. down-sample)
	// and do a final inverse FT
	if (ref_dim == 2)
		vol_out.resize(ori_size, ori_size);
	else
		vol_out.resize(ori_size, ori_size, ori_size);
	RCTOC(ReconTimer,ReconS_9);
	RCTIC(ReconTimer,ReconS_10);
	FourierTransformer transformer2;
	MultidimArray<Complex > Ftmp;
	transformer2.setReal(vol_out); // cannot use the first transformer because Fconv is inside there!!
	transformer2.getFourierAlias(Ftmp);
	RCTOC(ReconTimer,ReconS_10);
	RCTIC(ReconTimer,ReconS_11);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Ftmp)
	{
		if (kp * kp + ip * ip + jp * jp < r_max * r_max)
		{
			DIRECT_A3D_ELEM(Ftmp, k, i, j) = FFTW_ELEM(Fconv, kp * padding_factor, ip * padding_factor, jp * padding_factor);
		}
		else
		{
			DIRECT_A3D_ELEM(Ftmp, k, i, j) = 0.;
		}
	}
	RCTOC(ReconTimer,ReconS_11);
	RCTIC(ReconTimer,ReconS_13);
	CenterFFTbySign(Ftmp);
	RCTOC(ReconTimer,ReconS_13);
	RCTIC(ReconTimer,ReconS_12);
	// inverse FFT leaves result in vol_out
	transformer2.inverseFourierTransform();
	RCTOC(ReconTimer,ReconS_12);
	RCTIC(ReconTimer,ReconS_14);
	// Un-normalize FFTW (because original FFTs were done with the size of 2D FFTs)
	if (ref_dim==3)
		vol_out /= ori_size;
	RCTOC(ReconTimer,ReconS_14);
	RCTIC(ReconTimer,ReconS_15);
	// Mask out corners to prevent aliasing artefacts
	softMaskOutsideMap(vol_out);
	RCTOC(ReconTimer,ReconS_15);
	RCTIC(ReconTimer,ReconS_16);
	// Gridding correction for the blob
	RFLOAT normftblob = tab_ftblob(0.);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_out)
	{

		RFLOAT r = sqrt((RFLOAT)(k*k+i*i+j*j));
		RFLOAT rval = r / (ori_size * padding_factor);
		A3D_ELEM(vol_out, k, i, j) /= tab_ftblob(rval) / normftblob;
		//if (k==0 && i==0)
		//	std::cerr << " j= " << j << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << std::endl;
	}
	RCTOC(ReconTimer,ReconS_16);

#else

	// rather than doing the blob-convolution to downsample the data array, do a windowing operation:
	// This is the same as convolution with a SINC. It seems to give better maps.
	// Then just make the blob look as much as a SINC as possible....
	// The "standard" r1.9, m2 and a15 blob looks quite like a sinc until the first zero (perhaps that's why it is standard?)
	//for (RFLOAT r = 0.1; r < 10.; r+=0.01)
	//{
	//	RFLOAT sinc = sin(PI * r / padding_factor ) / ( PI * r / padding_factor);
	//	std::cout << " r= " << r << " sinc= " << sinc << " blob= " << blob_val(r, blob) << std::endl;
	//}

	// Now do inverse FFT and window to original size in real-space
	// Pass the transformer to prevent making and clearing a new one before clearing the one declared above....
	// The latter may give memory problems as detected by electric fence....
	RCTIC(ReconTimer,ReconS_17);
	windowToOridimRealSpace(transformer, vol_out, printTimes);
	RCTOC(ReconTimer,ReconS_17);

#endif

#ifdef DEBUG_RECONSTRUCT
	ttt()=vol_out;
	ttt.write("reconstruct_before_gridding_correction.spi");
#endif

	// Correct for the linear/nearest-neighbour interpolation that led to the data array
	RCTIC(ReconTimer,ReconS_18);

	griddingCorrect(vol_out);

	RCTOC(ReconTimer,ReconS_18);
	RCTIC(ReconTimer,ReconS_23);

	// Completely empty the transformer object
	transformer.cleanup();
	// Now can use extra mem to move data into smaller array space
	vol_out.shrinkToFit();

	RCTOC(ReconTimer,ReconS_23);

#ifdef TIMING
	if(printTimes)
		ReconTimer.printTimes(true);
#endif


#ifdef DEBUG_RECONSTRUCT
	std::cerr<<"done with reconstruct"<<std::endl;
#endif

	if (weight_out != 0)
	{
		weight_out->data = MultidimArray<RFLOAT>(1, ori_size, ori_size, ori_size/2+1);

		Image<RFLOAT> count(ori_size/2+1, ori_size, ori_size);
		count.data.initZeros();

		// downsample while considering padding:

		for (long int z = 0; z < Fweight.zdim; z++)
		for (long int y = 0; y < Fweight.ydim; y++)
		for (long int x = 0; x < Fweight.xdim; x++)
		{
			int xl = x;
			int yl = y < Fweight.ydim/2? y : y - Fweight.ydim;
			int zl = z < Fweight.zdim/2? z : z - Fweight.zdim;

			if (xl == Fweight.xdim - 1
			 || yl == Fweight.ydim/2 || yl == -Fweight.ydim/2 - 1
			 || zl == Fweight.zdim/2 || zl == -Fweight.zdim/2 - 1)
			{
				continue;
			}

			int xx = ROUND(xl / padding_factor);
			int yy = (ROUND(yl / padding_factor) + ori_size) % ori_size;
			int zz = (ROUND(zl / padding_factor) + ori_size) % ori_size;

			if (xx >= 0 && xx < ori_size/2+1
			 && yy >= 0 && yy < ori_size
			 && zz >= 0 && zz < ori_size)
			{
				DIRECT_A3D_ELEM(weight_out->data, zz, yy, xx) += DIRECT_A3D_ELEM(Fweight, z, y, x);
				DIRECT_A3D_ELEM(count.data, zz, yy, xx) += 1.0;
			}
		}

		const double pad3 = padding_factor * padding_factor * padding_factor;

		for (long int z = 0; z < ori_size; z++)
		for (long int y = 0; y < ori_size; y++)
		for (long int x = 0; x < ori_size/2 + 1; x++)
		{
			const RFLOAT c = DIRECT_A3D_ELEM(count.data, z, y, x);

			if (c > 0.0)
			{
				DIRECT_A3D_ELEM(weight_out->data, z, y, x) *= pad3/c;
			}
		}
	}
}

void BackProjector::reweightGrad() {

	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(data) {
		const int r2 = k * k + i * i + j * j;
		if (r2 < max_r2)
			A3D_ELEM(data, k, i, j) = A3D_ELEM(data, k, i, j) / XMIPP_MAX(1., A3D_ELEM(weight, k, i, j));
	}
}

void BackProjector::getFristMoment(
		MultidimArray<Complex> &mom,
		RFLOAT lambda)
{
	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);

	mom.setXmippOrigin();
	mom.xinit = 0;

	if (mom.sum() == 0.)
	{
		FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
		{
			const int r2 = k * k + i * i + j * j;
			if (r2 < max_r2)
				A3D_ELEM(mom, k, i, j) = A3D_ELEM(data, k, i, j);
		}
	}
	else
	{
		FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
		{
			const int r2 = k * k + i * i + j * j;
			if (r2 < max_r2)
			{
				Complex m = lambda * A3D_ELEM(mom, k, i, j) + (1 - lambda) * A3D_ELEM(data, k, i, j);
				A3D_ELEM(mom, k, i, j) = m;
			}
		}
	}
}

void BackProjector::getSecondMoment(
		MultidimArray<Complex> &mom,
		MultidimArray<Complex> &data_other,
		RFLOAT lambda)
{
	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);

	mom.setXmippOrigin();
	mom.xinit = 0;

	FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
	{
		const int r2 = k * k + i * i + j * j;
		if (r2 < max_r2)
		{
			Complex n = A3D_ELEM(data_other, k, i, j) - A3D_ELEM(data, k, i, j);
			Complex s = (A3D_ELEM(data_other, k, i, j) + A3D_ELEM(data, k, i, j)) / 2;

			A3D_ELEM(mom, k, i, j).real = lambda * A3D_ELEM(mom, k, i, j).real +
			                              (1 - lambda) * (norm(n)/(norm(s) + 1e-12));
			A3D_ELEM(mom, k, i, j).imag = 0.;
		}
	}
}

void BackProjector::applyMomenta(
		MultidimArray<Complex> &mom1_half1,
		MultidimArray<Complex> &mom1_half2,
		MultidimArray<Complex> &mom2)
{
	bool do_half = mom1_half2.nzyxdim == mom1_half1.nzyxdim;
	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);

	// Momentum application --------------------------------------------------------------------

	FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
	{
		const int r2 = k * k + i * i + j * j;
		if (r2 < max_r2)
		{
			// First moment
			Complex g = A3D_ELEM(mom1_half1, k, i, j);
			if (do_half)
				g = (g + A3D_ELEM(mom1_half2, k, i, j)) / 2.;


			// Second moment
			g /= sqrt(A3D_ELEM(mom2, k, i, j).real) + 1e-12;

			A3D_ELEM(data, k, i, j) = g;
		}
	}

	// FSC estimation -----------------------------------------------------------------------
	if (do_half)
	{
		mom1_noise_power.initZeros(ori_size / 2 + 1);

		MultidimArray<RFLOAT> counter;
		counter.initZeros(ori_size / 2 + 1);
		FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
				{
					const int r2 = k * k + i * i + j * j;
					if (r2 < max_r2)
					{
						int ires = ROUND(sqrt((RFLOAT) r2) / padding_factor);
						Complex diff = A3D_ELEM(mom1_half1, k, i, j) - A3D_ELEM(mom1_half2, k, i, j);
						DIRECT_A1D_ELEM(mom1_noise_power, ires) += sqrt(norm(diff));
						DIRECT_A1D_ELEM(counter, ires) += 1;
					}
				}
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(counter)
		{
			if (DIRECT_A1D_ELEM(counter, i) > 0.)
				DIRECT_A1D_ELEM(mom1_noise_power, i) /= DIRECT_A1D_ELEM(counter, i);
		}
	}
}

void BackProjector::reconstructGrad(
		MultidimArray<RFLOAT> &vol_out,
		const MultidimArray<RFLOAT> &fsc_spectrum,
		RFLOAT grad_stepsize,
		RFLOAT tau2_fudge,
		RFLOAT min_resol_shell,
		bool use_fsc,
		bool printTimes)
{
	const int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	RFLOAT oversampling_correction = (ref_dim == 3) ? (padding_factor * padding_factor * padding_factor) : (padding_factor * padding_factor);

	// Calculate Fourier transform of the current reference in the same way as done in the E-step
	MultidimArray<RFLOAT> dummy;
	Projector PPref(ori_size, interpolator, padding_factor, r_min_nn, data_dim);
	PPref.computeFourierTransformMap(vol_out, dummy, r_max*2, 1, false); // false means no gridding correction

//#define DEBUG_NGD
#ifdef DEBUG_NGD
	std::cerr << " oversampling_correction= " << oversampling_correction << std::endl;
	std::cerr << " subset_size= " << subset_size << std::endl;
	std::cerr << " pad_size= " << pad_size << " PPref.pad_size= " << PPref.pad_size << std::endl;
	std::cerr << " r_max= " << r_max << " PPref.r_max= " << PPref.r_max << std::endl;
	std::cerr << "PPref.data "; PPref.data.printShape(std::cerr);
	std::cerr << "data "; data.printShape(std::cerr);
	std::cerr << "weight "; weight.printShape(std::cerr);
	std::cerr << "Fconv "; Fconv.printShape(std::cerr);
#endif

//#define WRITE_DIFF
#ifdef WRITE_DIFF
	MultidimArray<Complex> Fdiff_cen(PPref.data);
#endif

    //---------------------- FSC estimate ----------------------//

	MultidimArray<RFLOAT> fsc_estimate;

	if (!use_fsc) {
		MultidimArray<RFLOAT> prev_power;
		MultidimArray<RFLOAT> counter;
		prev_power.initZeros(ori_size / 2 + 1);
		counter.initZeros(ori_size / 2 + 1);
		fsc_estimate.initZeros(ori_size / 2 + 1);

		// First get average Delta^2 from the accumulated gradient in data and weight
		FOR_ALL_ELEMENTS_IN_ARRAY3D(data) {
			const int r2 = k * k + i * i + j * j;
			if (r2 < max_r2) {
				int ires = ROUND(sqrt((RFLOAT) r2) / padding_factor);
				DIRECT_A1D_ELEM(prev_power, ires) += sqrt(norm(A3D_ELEM(PPref.data, k, i, j)));
				DIRECT_A1D_ELEM(counter, ires) += 1;
			}
		}

		RFLOAT threshold(0);

		if (mom1_noise_power.nzyxdim > 0)
		{
			//Average power spectra and calculate FSC estimate
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(counter)
			{
				RFLOAT prev = 0.;
				RFLOAT snr = 0.;

				if (DIRECT_A1D_ELEM(counter, i) > 0.)
				{
					prev = DIRECT_A1D_ELEM(prev_power, i) / DIRECT_A1D_ELEM(counter, i);
					DIRECT_A1D_ELEM(prev_power, i) = prev;

					if (DIRECT_A1D_ELEM(mom1_noise_power, i) > 0.)
						snr = 2 * tau2_fudge * prev / DIRECT_A1D_ELEM(mom1_noise_power, i);

					RFLOAT f = snr / (1 + snr);
					DIRECT_A1D_ELEM(fsc_estimate, i) = XMIPP_MIN(XMIPP_MAX(f, DIRECT_A1D_ELEM(fsc_spectrum, i)), 1);
				}
			}
		}
		else
			fsc_estimate += 1;

		tau2_fudge = 1;
	}
	else
		fsc_estimate = fsc_spectrum;

    //---------------------- Make gradient update ----------------------//

	FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
	{
		const int r2 = k * k + i * i + j * j;
		if (r2 < max_r2)
		{
			int ires = ROUND(sqrt((RFLOAT)r2) / padding_factor);
			RFLOAT fsc = DIRECT_A1D_ELEM(fsc_estimate, ires);

			Complex Fgrad = fsc * A3D_ELEM(data, k, i, j) - (1. - fsc) / tau2_fudge * A3D_ELEM(PPref.data, k, i, j);
			A3D_ELEM(PPref.data, k, i, j) += grad_stepsize * Fgrad;
		}
	}
#ifdef WRITE_DIFF
	{
		FourierTransformer transformer;
		transformer.setReal(vol_out);
		MultidimArray<Complex>& Ftmp = transformer.getFourierReference();
		Projector::decenter(data, Ftmp, max_r2);
		Image<RFLOAT> Idiff;
		windowToOridimRealSpace(transformer, Idiff(), printTimes);
		Idiff.write("diff.mrc");
	}
#endif

	// Set Fconv to the right size
	if (ref_dim == 2) {
		vol_out.setDimensions(pad_size, pad_size, 1, 1);
	}
	else {
		// Too costly to actually allocate the space
		// Trick transformer with the right dimensions
		vol_out.setDimensions(pad_size, pad_size, pad_size, 1);
	}

	FourierTransformer transformer;
	transformer.setReal(vol_out); // Fake set real. 1. Allocate space for Fconv 2. calculate plans.
	MultidimArray<Complex>& Fconv = transformer.getFourierReference();

	// Now do inverse FFT and window to original size in real-space
	// Pass the transformer to prevent making and clearing a new one before clearing the one declared above....
	// The latter may give memory problems as detected by electric fence....
	Projector::decenter(PPref.data, Fconv, max_r2);
//	RCTIC(ReconTimer,ReconS_17);
	windowToOridimRealSpace(transformer, vol_out, printTimes);
//	RCTOC(ReconTimer,ReconS_17);

#ifdef DEBUG_NGD
	Image<RFLOAT> Itmp;
	Itmp()=vol_out;
	Itmp.write("debug.spi");
#endif
}

void BackProjector::symmetrise(int nr_helical_asu, RFLOAT helical_twist, RFLOAT helical_rise, int threads)
{
	// First make sure the input arrays are obeying Hermitian symmetry,
	// which is assumed in the rotation operators of both helical and point group symmetry
	enforceHermitianSymmetry();

	// Then apply helical and point group symmetry (order irrelevant?)
	applyHelicalSymmetry(nr_helical_asu, helical_twist, helical_rise);

	applyPointGroupSymmetry(threads);
}

void BackProjector::enforceHermitianSymmetry()
{
	for (int iz = STARTINGZ(data); iz <=FINISHINGZ(data); iz++)
	{
		// Make sure all points are only included once.
		int starty = (iz < 0) ? 0 : 1;
		for (int iy = starty; iy <= FINISHINGY(data); iy++)
		{
			// I just need to sum the two points, not divide by 2!
			Complex fsum = (A3D_ELEM(data, iz, iy, 0) + conj(A3D_ELEM(data, -iz, -iy, 0)));
			A3D_ELEM(data, iz, iy, 0) = fsum;
			A3D_ELEM(data, -iz, -iy, 0) = conj(fsum);
			RFLOAT sum = (A3D_ELEM(weight, iz, iy, 0) + A3D_ELEM(weight, -iz, -iy, 0));
			A3D_ELEM(weight, iz, iy, 0) = sum;
			A3D_ELEM(weight, -iz, -iy, 0) = sum;
		}
	}
}

void BackProjector::applyHelicalSymmetry(int nr_helical_asu, RFLOAT helical_twist, RFLOAT helical_rise)
{
	if ( (nr_helical_asu < 2) || (ref_dim != 3) )
		return;

	int rmax2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);

	Matrix2D<RFLOAT> R(4, 4); // A matrix from the list
	MultidimArray<RFLOAT> sum_weight;
	MultidimArray<Complex > sum_data;
	RFLOAT x, y, z, fx, fy, fz, xp, yp, zp, r2;
	bool is_neg_x;
	int x0, x1, y0, y1, z0, z1;
	Complex d000, d001, d010, d011, d100, d101, d110, d111;
	Complex dx00, dx01, dx10, dx11, dxy0, dxy1, ddd;
	RFLOAT dd000, dd001, dd010, dd011, dd100, dd101, dd110, dd111;
	RFLOAT ddx00, ddx01, ddx10, ddx11, ddxy0, ddxy1;

	// First symmetry operator (not stored in SL) is the identity matrix
	sum_weight = weight;
	sum_data = data;
	int h_min = -nr_helical_asu/2;
	int h_max = -h_min + nr_helical_asu%2;
	for (int hh = h_min; hh < h_max; hh++)
	{
		if (hh != 0) // h==0 is done before the for loop (where sum_data = data)
		{
			RFLOAT rot_ang = hh * (-helical_twist);
			rotation3DMatrix(rot_ang, 'Z', R);
			R.setSmallValuesToZero(); // TODO: invert rotation matrix?

			// Loop over all points in the output (i.e. rotated, or summed) array
			FOR_ALL_ELEMENTS_IN_ARRAY3D(sum_weight)
			{
				x = (RFLOAT)j; // STARTINGX(sum_weight) is zero!
				y = (RFLOAT)i;
				z = (RFLOAT)k;
				r2 = x*x + y*y + z*z;
				if (r2 <= rmax2)
				{
					// coords_output(x,y) = A * coords_input (xp,yp)
					xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
					yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
					zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);

					// Only asymmetric half is stored
					if (xp < 0)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
						is_neg_x = true;
					}
					else
					{
						is_neg_x = false;
					}

					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					x0 = FLOOR(xp);
					fx = xp - x0;
					x1 = x0 + 1;

					y0 = FLOOR(yp);
					fy = yp - y0;
					y0 -=  STARTINGY(data);
					y1 = y0 + 1;

					z0 = FLOOR(zp);
					fz = zp - z0;
					z0 -= STARTINGZ(data);
					z1 = z0 + 1;

#ifdef CHECK_SIZE
					if (x0 < 0 || y0 < 0 || z0 < 0 ||
						x1 < 0 || y1 < 0 || z1 < 0 ||
						x0 >= XSIZE(data) || y0  >= YSIZE(data) || z0 >= ZSIZE(data) ||
						x1 >= XSIZE(data) || y1  >= YSIZE(data)  || z1 >= ZSIZE(data) 	)
					{
						std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
						std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
						data.printShape();
						REPORT_ERROR("BackProjector::applyPointGroupSymmetry: checksize!!!");
					}
#endif
					// First interpolate (complex) data
					d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
					d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
					d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
					d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
					d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
					d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
					d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
					d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

					dx00 = LIN_INTERP(fx, d000, d001);
					dx01 = LIN_INTERP(fx, d100, d101);
					dx10 = LIN_INTERP(fx, d010, d011);
					dx11 = LIN_INTERP(fx, d110, d111);
					dxy0 = LIN_INTERP(fy, dx00, dx10);
					dxy1 = LIN_INTERP(fy, dx01, dx11);

					// Take complex conjugated for half with negative x
					ddd = LIN_INTERP(fz, dxy0, dxy1);

					if (is_neg_x)
						ddd = conj(ddd);

					// Also apply a phase shift for helical translation along Z
					if (ABS(helical_rise) > 0.)
					{
						RFLOAT zshift = hh * helical_rise;
						zshift /= - ori_size * (RFLOAT)padding_factor;
						RFLOAT dotp = 2 * PI * (z * zshift);
						RFLOAT a = cos(dotp);
						RFLOAT b = sin(dotp);
						RFLOAT c = ddd.real;
						RFLOAT d = ddd.imag;
						RFLOAT ac = a * c;
						RFLOAT bd = b * d;
						RFLOAT ab_cd = (a + b) * (c + d);
						ddd = Complex(ac - bd, ab_cd - ac - bd);
					}
					// Accumulated sum of the data term
					A3D_ELEM(sum_data, k, i, j) += ddd;

					// Then interpolate (real) weight
					dd000 = DIRECT_A3D_ELEM(weight, z0, y0, x0);
					dd001 = DIRECT_A3D_ELEM(weight, z0, y0, x1);
					dd010 = DIRECT_A3D_ELEM(weight, z0, y1, x0);
					dd011 = DIRECT_A3D_ELEM(weight, z0, y1, x1);
					dd100 = DIRECT_A3D_ELEM(weight, z1, y0, x0);
					dd101 = DIRECT_A3D_ELEM(weight, z1, y0, x1);
					dd110 = DIRECT_A3D_ELEM(weight, z1, y1, x0);
					dd111 = DIRECT_A3D_ELEM(weight, z1, y1, x1);

					ddx00 = LIN_INTERP(fx, dd000, dd001);
					ddx01 = LIN_INTERP(fx, dd100, dd101);
					ddx10 = LIN_INTERP(fx, dd010, dd011);
					ddx11 = LIN_INTERP(fx, dd110, dd111);
					ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
					ddxy1 = LIN_INTERP(fy, ddx01, ddx11);

					A3D_ELEM(sum_weight, k, i, j) +=  LIN_INTERP(fz, ddxy0, ddxy1);

				} // end if r2 <= rmax2
			 } // end loop over all elements of sum_weight
		} // end if hh!=0
	} // end loop over hh

	data = sum_data;
	weight = sum_weight;
}

void BackProjector::applyPointGroupSymmetry(int threads)
{

//#define DEBUG_SYMM
#ifdef DEBUG_SYMM
	std::cerr << " SL.SymsNo()= " << SL.SymsNo() << std::endl;
	std::cerr << " SL.true_symNo= " << SL.true_symNo << std::endl;
#endif

	int rmax2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	if (SL.SymsNo() > 0 && ref_dim == 3)
	{
		Matrix2D<RFLOAT> L(4, 4), R(4, 4); // A matrix from the list
		MultidimArray<RFLOAT> sum_weight;
		MultidimArray<Complex > sum_data;

		// First symmetry operator (not stored in SL) is the identity matrix
		sum_weight = weight;
		sum_data = data;
		// Loop over all other symmetry operators
		for (int isym = 0; isym < SL.SymsNo(); isym++)
		{
			SL.get_matrices(isym, L, R);
#ifdef DEBUG_SYMM
			std::cerr << " isym= " << isym << " R= " << R << std::endl;
#endif
			// Loop over all points in the output (i.e. rotated, or summed) array

			#pragma omp parallel for num_threads(threads)
			for (long int k=STARTINGZ(sum_weight); k<=FINISHINGZ(sum_weight); k++)
			for (long int i=STARTINGY(sum_weight); i<=FINISHINGY(sum_weight); i++)
			for (long int j=STARTINGX(sum_weight); j<=FINISHINGX(sum_weight); j++)
			{
	        	RFLOAT x = (RFLOAT)j; // STARTINGX(sum_weight) is zero!
	        	RFLOAT y = (RFLOAT)i;
	        	RFLOAT z = (RFLOAT)k;
	        	RFLOAT r2 = x*x + y*y + z*z;

	        	if (r2 <= rmax2)
	        	{
	        		// coords_output(x,y) = A * coords_input (xp,yp)
					RFLOAT xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
					RFLOAT yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
					RFLOAT zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);

					bool is_neg_x;

					// Only asymmetric half is stored
					if (xp < 0)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
						is_neg_x = true;
					}
					else
					{
						is_neg_x = false;
					}

					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					int x0 = FLOOR(xp);
					RFLOAT fx = xp - x0;
					int x1 = x0 + 1;

					int y0 = FLOOR(yp);
					RFLOAT fy = yp - y0;
					y0 -=  STARTINGY(data);
					int y1 = y0 + 1;

					int z0 = FLOOR(zp);
					RFLOAT fz = zp - z0;
					z0 -= STARTINGZ(data);
					int z1 = z0 + 1;

#ifdef CHECK_SIZE
					if (x0 < 0 || y0 < 0 || z0 < 0 ||
						x1 < 0 || y1 < 0 || z1 < 0 ||
						x0 >= XSIZE(data) || y0  >= YSIZE(data) || z0 >= ZSIZE(data) ||
						x1 >= XSIZE(data) || y1  >= YSIZE(data)  || z1 >= ZSIZE(data) 	)
					{
						std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
						std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
						data.printShape();
						REPORT_ERROR("BackProjector::applyPointGroupSymmetry: checksize!!!");
					}
#endif
					// First interpolate (complex) data
					Complex d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
					Complex d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
					Complex d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
					Complex d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
					Complex d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
					Complex d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
					Complex d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
					Complex d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

					Complex dx00 = LIN_INTERP(fx, d000, d001);
					Complex dx01 = LIN_INTERP(fx, d100, d101);
					Complex dx10 = LIN_INTERP(fx, d010, d011);
					Complex dx11 = LIN_INTERP(fx, d110, d111);

					Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
					Complex dxy1 = LIN_INTERP(fy, dx01, dx11);

					// Take complex conjugated for half with negative x
					if (is_neg_x)
					{
						A3D_ELEM(sum_data, k, i, j) += conj(LIN_INTERP(fz, dxy0, dxy1));
					}
					else
					{
						A3D_ELEM(sum_data, k, i, j) += LIN_INTERP(fz, dxy0, dxy1);
					}

					// Then interpolate (real) weight
					RFLOAT dd000 = DIRECT_A3D_ELEM(weight, z0, y0, x0);
					RFLOAT dd001 = DIRECT_A3D_ELEM(weight, z0, y0, x1);
					RFLOAT dd010 = DIRECT_A3D_ELEM(weight, z0, y1, x0);
					RFLOAT dd011 = DIRECT_A3D_ELEM(weight, z0, y1, x1);
					RFLOAT dd100 = DIRECT_A3D_ELEM(weight, z1, y0, x0);
					RFLOAT dd101 = DIRECT_A3D_ELEM(weight, z1, y0, x1);
					RFLOAT dd110 = DIRECT_A3D_ELEM(weight, z1, y1, x0);
					RFLOAT dd111 = DIRECT_A3D_ELEM(weight, z1, y1, x1);

					RFLOAT ddx00 = LIN_INTERP(fx, dd000, dd001);
					RFLOAT ddx01 = LIN_INTERP(fx, dd100, dd101);
					RFLOAT ddx10 = LIN_INTERP(fx, dd010, dd011);
					RFLOAT ddx11 = LIN_INTERP(fx, dd110, dd111);

					RFLOAT ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
					RFLOAT ddxy1 = LIN_INTERP(fy, ddx01, ddx11);

					A3D_ELEM(sum_weight, k, i, j) += LIN_INTERP(fz, ddxy0, ddxy1);

	        	} // end if r2 <= rmax2

	        } // end loop over all elements of sum_weight

	    } // end loop over symmetry operators

	    data = sum_data;
	    weight = sum_weight;
	    // Average
	    // The division should only be done if we would search all (C1) directions, not if we restrict the angular search!
	    /*
	    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
	    {
	    	DIRECT_MULTIDIM_ELEM(data, n) = DIRECT_MULTIDIM_ELEM(sum_data, n) / (RFLOAT)(SL.SymsNo() + 1);
	    	DIRECT_MULTIDIM_ELEM(weight, n) = DIRECT_MULTIDIM_ELEM(sum_weight, n) / (RFLOAT)(SL.SymsNo() + 1);
	    }
	    */
	}

}

void BackProjector::convoluteBlobRealSpace(FourierTransformer &transformer, bool do_mask)
{

	MultidimArray<RFLOAT> Mconv;
	int padhdim = pad_size / 2;

	// Set up right dimension of real-space array
	// TODO: resize this according to r_max!!!
	if (ref_dim==2)
		Mconv.reshape(pad_size, pad_size);
	else
		Mconv.reshape(pad_size, pad_size, pad_size);

	// inverse FFT
	transformer.setReal(Mconv);
	transformer.inverseFourierTransform();

	// Blob normalisation in Fourier space
	RFLOAT normftblob = tab_ftblob(0.);

	// TMP DEBUGGING
	//struct blobtype blob;
	//blob.order = 0;
	//blob.radius = 1.9 * padding_factor;
	//blob.alpha = 15;

    // Multiply with FT of the blob kernel
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Mconv)
    {
		int kp = (k < padhdim) ? k : k - pad_size;
		int ip = (i < padhdim) ? i : i - pad_size;
		int jp = (j < padhdim) ? j : j - pad_size;
    	RFLOAT rval = sqrt ( (RFLOAT)(kp * kp + ip * ip + jp * jp) ) / (ori_size * padding_factor);
    	//if (kp==0 && ip==0 && jp > 0)
		//	std::cerr << " jp= " << jp << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << " ori_size/2= " << ori_size/2 << std::endl;
    	// In the final reconstruction: mask the real-space map beyond its original size to prevent aliasing ghosts
    	// Note that rval goes until 1/2 in the oversampled map
    	if (do_mask && rval > 1./(2. * padding_factor))
    		DIRECT_A3D_ELEM(Mconv, k, i, j) = 0.;
    	else
    		DIRECT_A3D_ELEM(Mconv, k, i, j) *= (tab_ftblob(rval) / normftblob);
    }

    // forward FFT to go back to Fourier-space
    transformer.FourierTransform();
}

void BackProjector::windowToOridimRealSpace(FourierTransformer &transformer, MultidimArray<RFLOAT> &Mout, bool printTimes)
{

#ifdef TIMING
	Timer OriDimTimer;
	int OriDim1  = OriDimTimer.setNew(" OrD1_getFourier ");
	int OriDim2  = OriDimTimer.setNew(" OrD2_windowFFT ");
	int OriDim3  = OriDimTimer.setNew(" OrD3_reshape ");
	int OriDim4  = OriDimTimer.setNew(" OrD4_setReal ");
	int OriDim5  = OriDimTimer.setNew(" OrD5_invFFT ");
	int OriDim6  = OriDimTimer.setNew(" OrD6_centerFFT ");
	int OriDim7  = OriDimTimer.setNew(" OrD7_window ");
	int OriDim8  = OriDimTimer.setNew(" OrD8_norm ");
	int OriDim9  = OriDimTimer.setNew(" OrD9_softMask ");
#endif

	RCTIC(OriDimTimer,OriDim1);
	MultidimArray<Complex>& Fin = transformer.getFourierReference();
	RCTOC(OriDimTimer,OriDim1);
	RCTIC(OriDimTimer,OriDim2);
	MultidimArray<Complex > Ftmp;
	// Size of padded real-space volume
	int padoridim = ROUND(padding_factor * ori_size);
	// make sure padoridim is even
	padoridim += padoridim%2;
	RFLOAT normfft;

//#define DEBUG_WINDOWORIDIMREALSPACE
#ifdef DEBUG_WINDOWORIDIMREALSPACE
	Image<RFLOAT> tt;
	tt().reshape(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin)
	{
		DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
	}
	tt.write("windoworidim_Fin.spi");
#endif

    // Resize incoming complex array to the correct size
	windowFourierTransform(Fin, padoridim);
	RCTOC(OriDimTimer,OriDim2);
	RCTIC(OriDimTimer,OriDim3);
 	if (ref_dim == 2)
	{
		Mout.reshape(padoridim, padoridim);
		if (data_dim == 2)
			normfft = (RFLOAT)(padding_factor * padding_factor);
		else
			normfft = (RFLOAT)(padding_factor * padding_factor * ori_size);
	}
	else
	{
		Mout.reshape(padoridim, padoridim, padoridim);
		if (data_dim == 3)
			normfft = (RFLOAT)(padding_factor * padding_factor * padding_factor);
		else
			normfft = (RFLOAT)(padding_factor * padding_factor * padding_factor * ori_size);
	}
	Mout.setXmippOrigin();
	RCTOC(OriDimTimer,OriDim3);

#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt().reshape(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin)
	{
		DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
	}
	tt.write("windoworidim_Fresized.spi");
#endif

	// Shift the map back to its origin
	RCTIC(OriDimTimer,OriDim6);
	CenterFFTbySign(Fin);
	RCTOC(OriDimTimer,OriDim6);

	// Do the inverse FFT
	RCTIC(OriDimTimer,OriDim4);
	transformer.setReal(Mout);
	RCTOC(OriDimTimer,OriDim4);
	RCTIC(OriDimTimer,OriDim5);
#ifdef TIMING
	if(printTimes)
		std::cout << std::endl << "FFTrealDims = (" << transformer.fReal->xdim << " , " << transformer.fReal->ydim << " , " << transformer.fReal->zdim << " ) " << std::endl;
#endif
	transformer.inverseFourierTransform();
	RCTOC(OriDimTimer,OriDim5);
	//transformer.inverseFourierTransform(Fin, Mout);
	Fin.clear();
	transformer.fReal = NULL; // Make sure to re-calculate fftw plan
	Mout.setXmippOrigin();

#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt()=Mout;
	tt.write("windoworidim_Munwindowed.spi");
#endif

	// Window in real-space
	RCTIC(OriDimTimer,OriDim7);
	if (ref_dim==2)
	{
		Mout.window(FIRST_XMIPP_INDEX(ori_size), FIRST_XMIPP_INDEX(ori_size),
				       LAST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size));
	}
	else
	{
		Mout.window(FIRST_XMIPP_INDEX(ori_size), FIRST_XMIPP_INDEX(ori_size), FIRST_XMIPP_INDEX(ori_size),
				       LAST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size));
	}
	Mout.setXmippOrigin();
	RCTOC(OriDimTimer,OriDim7);
	// Normalisation factor of FFTW
	// The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
	RCTIC(OriDimTimer,OriDim8);
	Mout /= normfft;
	RCTOC(OriDimTimer,OriDim8);
#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt()=Mout;
	tt.write("windoworidim_Mwindowed.spi");
#endif

	// Mask out corners to prevent aliasing artefacts
	RCTIC(OriDimTimer,OriDim9);
	softMaskOutsideMap(Mout);
	RCTOC(OriDimTimer,OriDim9);

#ifdef DEBUG_WINDOWORIDIMREALSPACE
	tt()=Mout;
	tt.write("windoworidim_Mwindowed_masked.spi");
	FourierTransformer ttf;
	ttf.FourierTransform(Mout, Fin);
	tt().resize(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin)
	{
		DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
	}
	tt.write("windoworidim_Fnew.spi");
#endif

#ifdef TIMING
    if(printTimes)
    	OriDimTimer.printTimes(true);
#endif
}
