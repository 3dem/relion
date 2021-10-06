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
/***************************************************************************
 *
 * Authors: Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 e You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef _CTF_HH
#define _CTF_HH

#include <src/multidim_array.h>
#include <src/metadata_table.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <map>


class CTF
{
protected:

	// Different constants
	RFLOAT K1;
	RFLOAT K2;
	RFLOAT K3;
	RFLOAT K4;
	RFLOAT K5;

	// Astigmatism stored in symmetrical matrix form
	RFLOAT Axx, Axy, Ayy;

	// Azimuthal angle in radians
	RFLOAT rad_azimuth;

	// defocus_average = (defocus_u + defocus_v)/2
	RFLOAT defocus_average;

	// defocus_deviation = (defocus_u - defocus_v)/2
	RFLOAT defocus_deviation;

	// Pointer to observation model kept after a call to readByGroup() to enable
	// caching of symmetric aberrations (CTF instances can be reallocated for each particle,
	// while the same obs. model lives for the entire duration of the program)
	ObservationModel* obsModel;
	int opticsGroup;

public:

	/// Accelerating Voltage (in KiloVolts)
	RFLOAT kV;

	/// Defocus in U (in Angstroms). Positive values are underfocused
	RFLOAT DeltafU;

	/// Defocus in V (in Angstroms). Postive values are underfocused
	RFLOAT DeltafV;

	/// Azimuthal angle (between X and U) in degrees
	RFLOAT azimuthal_angle;

	// Electron wavelength (Amstrongs)
	RFLOAT lambda;

	// Radius of the aperture (in micras)
	// RFLOAT aperture;
	// Spherical aberration (in milimeters). Typical value 5.6
	RFLOAT Cs;

	/// Chromatic aberration (in milimeters). Typical value 2
	RFLOAT Ca;

	/** Mean energy loss (eV) due to interaction with sample.
	    Typical value 1*/
	RFLOAT espr;

	/// Objective lens stability (deltaI/I) (ppm). Typical value 1
	RFLOAT ispr;

	/// Convergence cone semiangle (in mrad). Typical value 0.5
	RFLOAT alpha;

	/// Longitudinal mechanical displacement (Angstrom). Typical value 100
	RFLOAT DeltaF;

	/// Transversal mechanical displacement (Angstrom). Typical value 3
	RFLOAT DeltaR;

	/// Amplitude contrast. Typical values 0.07 for cryo, 0.2 for negative stain
	RFLOAT Q0;

	// B-factor fall-off
	RFLOAT Bfac;

	// Overall scale-factor of CTF
	RFLOAT scale;

	// Phase-shift from a phase-plate (in rad)
	RFLOAT phase_shift;

	/** Empty constructor. */
	CTF() :
		kV(200), DeltafU(0), DeltafV(0), azimuthal_angle(0), phase_shift(0),
		Cs(0), Bfac(0), Q0(0), scale(1), obsModel(0), opticsGroup(0)
	{}

	// Read CTF parameters from particle table partMdt and optics table opticsMdt.
	void readByGroup(const MetaDataTable &partMdt, ObservationModel* obs, long int particle = -1);

	void readValue(EMDLabel label, RFLOAT& dest, RFLOAT defaultVal, long int particle, int opticsGroup,
	               const MetaDataTable& partMdt, const ObservationModel* obs);

	/** Read CTF parameters from MetaDataTables MD1 and MD2 (deprecated)
	 * If a parameter is not found in MD1 it is tried to be read from MD2
	 * Only if it is also not present in the second then a default value is used
	 * This is useful if micrograph-specific parameters are stored in a separate MD from the image-specific parameters
	 */
	void read(const MetaDataTable &MD1, const MetaDataTable &MD2, long int objectID = -1);

	/** Just set all values explicitly */
	void setValues(RFLOAT _defU, RFLOAT _defV, RFLOAT _defAng,
	               RFLOAT _voltage, RFLOAT _Cs, RFLOAT _Q0, RFLOAT _Bfac, RFLOAT _scale = 1., RFLOAT _phase_shift = 0.);

	/** Set all values explicitly in 3.1 */
	void setValuesByGroup(ObservationModel* obs, int opticsGroup,
	                      RFLOAT _defU, RFLOAT _defV, RFLOAT _defAng,
	                      RFLOAT _Bfac = 0.0, RFLOAT _scale = 1.0, RFLOAT _phase_shift = 0.0);


	/** Read from a single MetaDataTable */
	void read(const MetaDataTable &MD);

	/** Write to MetaDataTable. */
	void write(MetaDataTable &MD);

	/** Write to output. */
	void write(std::ostream &out);

	/// Set up the CTF object, read parameters from MetaDataTables with micrograph and particle information
	/// If no MDmic is provided or it does not contain certain parameters, these parameters are tried to be read from MDimg
	void initialise();

	/// Compute CTF at (U,V). Continuous frequencies
	inline RFLOAT getCTF(RFLOAT X, RFLOAT Y,
	                     bool do_abs = false, bool do_only_flip_phases = false,
	                     bool do_intact_until_first_peak = false, bool do_damping = true,
	                     double gammaOffset = 0.0, bool do_intact_after_first_peak = false) const
	{
		if (obsModel != 0 && obsModel->hasMagMatrices)
		{
			const Matrix2D<RFLOAT>& M = obsModel->getMagMatrix(opticsGroup);
			RFLOAT Xd = M(0,0) * X + M(0,1) * Y;
			RFLOAT Yd = M(1,0) * X + M(1,1) * Y;

			X = Xd;
			Y = Yd;
		}

		RFLOAT u2 = X * X + Y * Y;
		RFLOAT u4 = u2 * u2;

		// if (u2>=ua2) return 0;
		//RFLOAT deltaf = getDeltaF(X, Y);
		//RFLOAT gamma = K1 * deltaf * u2 + K2 * u4 - K5 - K3 + gammaOffset;
		RFLOAT gamma = K1 * (Axx*X*X + 2.0*Axy*X*Y + Ayy*Y*Y) + K2 * u4 - K5 - K3 + gammaOffset;

		RFLOAT retval;

		if ((do_intact_until_first_peak && ABS(gamma) < PI/2.) ||
		    (do_intact_after_first_peak && ABS(gamma) > PI/2.))
		{
			retval = 1.;
		}
		else
		{
			retval = -sin(gamma);
		}

		if (do_damping)
		{
			RFLOAT E = exp(K4 * u2); // B-factor decay (K4 = -Bfac/4);
			retval *= E;
		}

		if (do_abs)
		{
			retval = ABS(retval);
		}
		else if (do_only_flip_phases)
		{
			retval = (retval < 0.) ? -1. : 1.;
		}

		retval *= scale;

		// SHWS 25-2-2019: testing a new idea to improve code stability
		// Don't allow very small values of CTF to prevent division by zero in GPU code
		if (fabs(retval) < 1e-8)
		{
			retval = SGN(retval) * 1e-8;
		}

		return retval;
	}

	RFLOAT getLowOrderGamma(RFLOAT X, RFLOAT Y) const;

	// compute the local frequency of the ctf
	// (i.e. the radial slope of 'double gamma' in getCTF())
	// -- deprecated, use getGammaGrad().length()
	RFLOAT getCtfFreq(RFLOAT X, RFLOAT Y);

	gravis::t2Vector<RFLOAT> getGammaGrad(RFLOAT X, RFLOAT Y) const;

	inline Complex getCTFP(RFLOAT X, RFLOAT Y, bool is_positive, double gammaOffset = 0.0) const
	{
		if (obsModel != 0 && obsModel->hasMagMatrices)
		{
			const Matrix2D<RFLOAT>& M = obsModel->getMagMatrix(opticsGroup);
			RFLOAT Xd = M(0,0) * X + M(0,1) * Y;
			RFLOAT Yd = M(1,0) * X + M(1,1) * Y;

			X = Xd;
			Y = Yd;
		}

		RFLOAT u2 = X * X + Y * Y;
		RFLOAT u4 = u2 * u2;

		RFLOAT gamma = K1 * (Axx*X*X + 2.0*Axy*X*Y + Ayy*Y*Y) + K2 * u4 - K5 - K3 + gammaOffset;

		RFLOAT sinx, cosx;
#ifdef RELION_SINGLE_PRECISION
		SINCOSF( gamma, &sinx, &cosx );
#else
		SINCOS( gamma, &sinx, &cosx );
#endif
		Complex retval;
		retval.real = -sinx;
		retval.imag = (is_positive) ? cosx : -cosx;

		return retval;
	}

	/// Compute Deltaf at a given direction (no longer used by getCTF)
	inline RFLOAT getDeltaF(RFLOAT X, RFLOAT Y) const
	{
		if (ABS(X) < XMIPP_EQUAL_ACCURACY &&
		    ABS(Y) < XMIPP_EQUAL_ACCURACY)
			return 0;

		RFLOAT ellipsoid_ang = atan2(Y, X) - rad_azimuth;
		/*
		* For a derivation of this formulae confer
		* Principles of Electron Optics page 1380
		* in particular term defocus and twofold axial astigmatism
		* take into account that a1 and a2 are the coefficient
		* of the zernike polynomials difference of defocus at 0
		* and at 45 degrees. In this case a2=0
		*/
		RFLOAT cos_ellipsoid_ang_2 = cos(2*ellipsoid_ang);
		return (defocus_average + defocus_deviation*cos_ellipsoid_ang_2);

	}

	/// Generate (Fourier-space, i.e. FFTW format) image with all CTF values.
	/// The dimensions of the result array should have been set correctly already
	//  FFTW format means the Nyquist component is on the positive side.
	//   e.g. for N = 6,  kx = [0, 1, 2, 3], ky = [0, 1, 2, 3, -2, -1]
	//  Unfortunately, codes in jaz/ use different convention.
	//   e.g. for N - 6, kx = [0, 1, 2, 3], ky = [0, 1, 2, -3, -2, -1]
	//  TODO: FIXME: Thus, the returned values at Nyquist are WRONG!!!
	void getFftwImage(MultidimArray < RFLOAT > &result, int orixdim, int oriydim, RFLOAT angpix,
	                  bool do_abs = false, bool do_only_flip_phases = false, bool do_intact_until_first_peak = false,
	                  bool do_damping = true, bool do_ctf_padding = false, bool do_intact_after_first_peak = false) const;

	// Get a complex image with the CTFP/Q values, where the angle is in degrees between the Y-axis and the CTFP/Q sector line
	void getCTFPImage(MultidimArray<Complex> &result, int orixdim, int oriydim, RFLOAT angpix,
	                  bool is_positive, float angle);

	/// Generate a centered image (with hermitian symmetry)
	/// The dimensions of the result array should have been set correctly already
	void getCenteredImage(MultidimArray < RFLOAT > &result, RFLOAT angpix,
	                      bool do_abs = false, bool do_only_flip_phases = false, bool do_intact_until_first_peak = false,
	                      bool do_damping = true, bool do_intact_after_first_peak = false);

	/// Generate a 1D profile along defocusAngle
	/// The dimensions of the result array should have been set correctly already, i.e. at the image size!
	void get1DProfile(MultidimArray < RFLOAT > &result, RFLOAT angle, RFLOAT angpix,
	                  bool do_abs = false, bool do_only_flip_phases = false, bool do_intact_until_first_peak = false, 
	                  bool do_damping = true, bool do_intact_after_first_peak = false);

	// Calculate weight W for Ewald-sphere curvature correction: apply this to the result from getFftwImage
	void applyWeightEwaldSphereCurvature(MultidimArray<RFLOAT>& result, int orixdim, int oriydim,
	                                     RFLOAT angpix, RFLOAT particle_diameter);

	void applyWeightEwaldSphereCurvature_new(MultidimArray<RFLOAT>& result, int orixdim, int oriydim,
	                                         RFLOAT angpix, RFLOAT particle_diameter);

	// Calculate weight W for Ewald-sphere curvature correction: apply this to the result from getFftwImage
	void applyWeightEwaldSphereCurvature_noAniso(MultidimArray < RFLOAT > &result, int orixdim, int oriydim, RFLOAT angpix, RFLOAT particle_diameter);

	void applyEwaldMask(RawImage<RFLOAT>& result, int orixdim, int oriydim, RFLOAT angpix, RFLOAT particle_diameter);

	
	std::vector<double> getK() const;
	double getAxx() const;
	double getAxy() const;
	double getAyy() const;
	
	
	
	// Methods using the new data types from 2019/2020 in src/jaz/
	
	template <typename T>
	void draw(
		int w0, int h0, double angpix,
		const BufferedImage<double>* gammaOffset, T* dest) const
	{
		const int wh = w0 / 2 + 1;
		
		double xs = w0 * angpix;
		double ys = h0 * angpix;
		
		if (gammaOffset)
		{
			if (gammaOffset->ydim != h0)
			{
				REPORT_ERROR_STR(
					"CTF::draw: wrong cached gamma-offset size. Box size: "
					<< w0 << ", cache size: " << gammaOffset->ydim);
			}

			for (int y = 0; y < h0; y++)
			for (int x = 0; x < wh; x++)
			{
				double xx = x / xs;
				double yy = (y < h0/2? y : y - h0) / ys;

				dest[y*wh + x] = getCTF(xx, yy, false, false, false, true, (*gammaOffset)(x,y));
			}
		}
		else
		{
			for (int y = 0; y < h0; y++)
			for (int x = 0; x < wh; x++)
			{
				double xx = x / xs;
				double yy = (y < h0/2? y : y - h0) / ys;

				dest[y*wh + x] = getCTF(xx,yy, false, false, false, true);
			}
		}
	}
	
	template <typename T>
	void draw_fast(int w0, int h0, double angpix,
				   const BufferedImage<double>* gammaOffset, T* dest) const
	{
		const int wh = w0 / 2 + 1;
		
		const double xs = w0 * angpix;
		const double ys = h0 * angpix;
		
		const double r2_max = wh * wh / (xs * xs);
		
		if (gammaOffset)
		{
			if (gammaOffset->ydim != h0)
			{
				REPORT_ERROR_STR(
					"CTF::draw_fast: wrong cached gamma-offset size. Box size: "
					<< w0 << ", cache size: " << gammaOffset->ydim);
			}

			for (int y = 0; y < h0; y++)
			for (int x = 0; x < wh; x++)
			{
				const double xx = x / xs;
				const double yy = (y < h0/2? y : y - h0) / ys;

				const double u2 = xx * xx + yy * yy;

				if (u2 < r2_max)
				{
					double u4 = u2 * u2;

					const double gamma = K1 * (Axx*xx*xx + 2.0*Axy*xx*yy + Ayy*yy*yy) + K2 * u4 - K5 - K3 + (*gammaOffset)(x,y);

					dest[y*wh + x] = -scale*sin(gamma);
				}
			}
		}
		else
		{
			for (int y = 0; y < h0; y++)
			for (int x = 0; x < wh; x++)
			{
				const double xx = x / xs;
				const double yy = (y < h0/2? y : y - h0) / ys;

				const double u2 = xx * xx + yy * yy;

				if (u2 < r2_max)
				{
					double u4 = u2 * u2;

					const double gamma = K1 * (Axx*xx*xx + 2.0*Axy*xx*yy + Ayy*yy*yy) + K2 * u4 - K5 - K3;

					dest[y*wh + x] = -sin(gamma);
				}
			}
		}
	}
	
	template <typename T>
	void drawGamma(int w0, int h0, double angpix, T* dest) const
	{
		const int wh = w0 / 2 + 1;
		
		double xs = w0 * angpix;
		double ys = h0 * angpix;
		
		if (obsModel == 0 || !obsModel->hasEvenZernike)
		{
			for (int y = 0; y < h0; y++)
			for (int x = 0; x < wh; x++)
			{
				double xx = x / xs;
				double yy = (y < h0/2? y : y - h0) / ys;
		
				dest[y*wh + x] = getLowOrderGamma(xx,yy);
			}
		}
		else
		{
			const BufferedImage<RFLOAT>& gammaOffset = obsModel->getGammaOffset(opticsGroup, h0);
			
			for (int y = 0; y < h0; y++)
			for (int x = 0; x < wh; x++)
			{
				double xx = x / xs;
				double yy = (y < h0/2? y : y - h0) / ys;
		
				dest[y*wh + x] = getLowOrderGamma(xx,yy) + gammaOffset(x,y);
			}
		}
	}
	
	void setDefocusMatrix(double axx, double axy, double ayy);
			
};
#endif
