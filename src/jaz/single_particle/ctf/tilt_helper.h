/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#ifndef TILT_REFINEMENT_H
#define TILT_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class TiltOptimization : public Optimization
{
    public:

        TiltOptimization(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double angpix,
		const Matrix2D<RFLOAT>& mag,
                bool L1 = false,
                bool anisotropic = false);

        double f(const std::vector<double>& x, void* tempStorage) const;

    private:

        const Image<Complex>& xy;
        const Image<RFLOAT>& weight;

        const double angpix;
	const Matrix2D<RFLOAT>& mag;
        const bool L1, anisotropic;
};

class BasisOptimisation : public Optimization
{
    public:

        BasisOptimisation(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                const std::vector<Image<RFLOAT>>& basis,
                bool L1 = false);

        double f(const std::vector<double>& x, void* tempStorage) const;

		void* allocateTempStorage() const;
		void deallocateTempStorage(void* ts) const;

    private:

		const int w, h, cc;
        const Image<Complex>& xy;
        const Image<RFLOAT>& weight;
		const std::vector<Image<RFLOAT>>& basis;
        const bool L1;
};

class AnisoBasisOptimisation : public Optimization
{
    public:

        AnisoBasisOptimisation(
                const Image<Complex>& xy,
				const Image<RFLOAT>& weight0,
				const Image<RFLOAT>& Axx,
				const Image<RFLOAT>& Axy,
				const Image<RFLOAT>& Ayy,
                const std::vector<Image<RFLOAT>>& basis,
                bool L1 = false);

        double f(const std::vector<double>& x, void* tempStorage) const;

		void* allocateTempStorage() const;
		void deallocateTempStorage(void* ts) const;

    private:

		const int w, h, cc;
        const Image<Complex>& xy;
        const Image<RFLOAT>& weight0, &Axx, &Axy, &Ayy;
		const std::vector<Image<RFLOAT>>& basis;
        const bool L1;
};

class TiltHelper
{
    public:

        static void updateTiltShift(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                CTF& ctf, double angpix,
                Image<Complex>& xyDest,
                Image<RFLOAT>& wDest,
				bool do_ctf_padding = false);

        static void updateTiltShiftPar(const Image<Complex>& prediction,
                const Image<Complex>& observation,
                CTF &ctf, double angpix,
                Image<Complex>& xyDest,
                Image<RFLOAT>& wDest,
				bool do_ctf_padding = false);

        static void fitTiltShift(
                const Image<RFLOAT>& phase,
                const Image<RFLOAT>& weight,
                double Cs, double lambda,
		double angpix, const Matrix2D<RFLOAT>& mag,
                double* shift_x, double* shift_y,
                double* tilt_x, double* tilt_y,
                Image<RFLOAT>* fit);

        static void optimizeTilt(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double Cs, double lambda,
		double angpix, const Matrix2D<RFLOAT>& mag,
                bool L1,
                double shift0_x, double shift0_y,
                double tilt0_x, double tilt0_y,
                double* shift_x, double* shift_y,
                double* tilt_x, double* tilt_y,
                Image<RFLOAT>* fit);


		static std::vector<double> fitOddZernike(
				const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
		double angpix, const Matrix2D<RFLOAT>& mag,
		int n_max,
		Image<RFLOAT>* fit = 0);

		static std::vector<double> optimiseOddZernike(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double angpix, const Matrix2D<RFLOAT>& mag,
		int n_max,
                const std::vector<double>& coeffs,
                Image<RFLOAT>* fit);

		// Nyquist X is positive, Y is negative (non-FFTW!!)
		static std::vector<Image<RFLOAT>> computeOddZernike(
				int s, double angpix, const Matrix2D<RFLOAT>& mag, int n_max);

		static Image<RFLOAT> plotOddZernike(
				const std::vector<double>& coeffs, int s,
				double angpix, const Matrix2D<RFLOAT>& mag);

		static Image<RFLOAT> plotTilt(
				double tx, double ty, int s,
				double angpix, const Matrix2D<RFLOAT>& mag,
				double Cs, double lambda);



		static std::vector<double> fitEvenZernike(
				const Image<RFLOAT>& phase,
		const Image<RFLOAT>& weight,
				double angpix, const Matrix2D<RFLOAT>& mag,
				int n_max, Image<RFLOAT>* fit = 0);

		static std::vector<double> optimiseEvenZernike(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double angpix, const Matrix2D<RFLOAT>& mag,
		int n_max,
                const std::vector<double>& coeffs,
                Image<RFLOAT>* fit);

		static std::vector<double> optimiseEvenZernike(
                const Image<Complex>& xy,
		const Image<RFLOAT>& weight0,
		const Image<RFLOAT>& Axx,
		const Image<RFLOAT>& Axy,
		const Image<RFLOAT>& Ayy,
                double angpix, const Matrix2D<RFLOAT>& mag,
		int n_max,
                const std::vector<double>& coeffs,
                Image<RFLOAT>* fit);

		// Nyquist X is positive, Y is negative (non-FFTW!!)
		static std::vector<Image<RFLOAT>> computeEvenZernike(
				int s,
				double angpix, const Matrix2D<RFLOAT>& mag,
				int n_max);

		static void extractTilt(
				std::vector<double>& oddZernikeCoeffs,
				double& tilt_x, double& tilt_y,
				double Cs, double lambda);

		static void insertTilt(
				std::vector<double>& oddZernikeCoeffs,
				double tilt_x, double tilt_y,
				double Cs, double lambda);


		static std::vector<double> fitBasisLin(
				const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
				const std::vector<Image<RFLOAT>>& basis);

		static std::vector<double> fitBasisLin(
				const Image<RFLOAT>& phase,
                const Image<RFLOAT>& weight,
				const std::vector<Image<RFLOAT>>& basis);

		static std::vector<double> optimiseBasis(
				const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
				const std::vector<Image<RFLOAT>>& basis,
				const std::vector<double>& initial);

		static std::vector<double> optimiseBasis(
				const Image<Complex>& xy,
				const Image<RFLOAT>& weight0,
				const Image<RFLOAT>& Axx,
				const Image<RFLOAT>& Axy,
				const Image<RFLOAT>& Ayy,
				const std::vector<Image<RFLOAT>>& basis,
				const std::vector<double>& initial);

        static void optimizeAnisoTilt(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double Cs, double lambda,
		double angpix, const Matrix2D<RFLOAT>& mag,
                bool L1,
                double shift0_x, double shift0_y,
                double tilt0_x, double tilt0_y,
                double* shift_x, double* shift_y,
                double* tilt_x, double* tilt_y,
                double* tilt_xx, double* tilt_xy, double* tilt_yy,
                Image<RFLOAT>* fit);

        static void drawPhaseShift(
                double shift_x, double shift_y,
                double tilt_x, double tilt_y,
                int w, int h, double as,
                const Matrix2D<RFLOAT>& mag,
                Image<RFLOAT>* tgt);

        static void drawPhaseShift(
                double shift_x, double shift_y,
                double tilt_x, double tilt_y,
                double tilt_xx, double tilt_xy, double tilt_yy,
                int w, int h, double as,
                const Matrix2D<RFLOAT>& mag,
                Image<RFLOAT>* tgt);

};

#endif
