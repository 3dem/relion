#ifndef ZERNIKE_HELPER_H
#define ZERNIKE_HELPER_H

#include "Zernike.h"
#include "tensor2x2.h"
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/optimization/optimization.h>

class CTF;

class ZernikeHelper
{
	public:
		
		
		static BufferedImage<double> computeEvenZernike(
				int s,
				double pixelSize, 
				const gravis::d2Matrix& mag,
				int n_max);
		
		static std::vector<double> fitEvenZernike(
				const BufferedImage<double>& phase,
				const BufferedImage<double>& weight,
				double pixelSize, const gravis::d2Matrix& mag,
				int n_max, BufferedImage<double>* fit = 0);
		
		static std::vector<double> optimiseEvenZernike(
				const BufferedImage<dComplex>& xy,
				const BufferedImage<Tensor2x2<double>>& A,
				double pixelSize, const gravis::d2Matrix& mag,
				int n_max,
				const std::vector<double>& coeffs,
				BufferedImage<double>* fit);
		
		
		
		static BufferedImage<double> computeOddZernike(
				int s,
				double pixelSize, 
				const gravis::d2Matrix& mag,
				int n_max);
		
		static std::vector<double> fitOddZernike(
				const BufferedImage<dComplex>& xy,
				const BufferedImage<double>& weight,
				double pixelSize, const gravis::d2Matrix& mag,
				int n_max,
				BufferedImage<double>* fit = 0);
		
		static std::vector<double> optimiseOddZernike(
				const BufferedImage<dComplex>& xy,
				const BufferedImage<double>& weight,
				double pixelSize, const gravis::d2Matrix& mag,
				int n_max,
				const std::vector<double>& coeffs,
				BufferedImage<double>* fit);
		
		
		
		static std::vector<double> fitBasisLin(
				const BufferedImage<dComplex>& xy,
				const BufferedImage<double>& weight,
				const BufferedImage<double>& basis);

		static std::vector<double> fitBasisLin(
				const BufferedImage<double>& phase,
				const BufferedImage<double>& weight,
				const BufferedImage<double>& basis);
		
		
		
		static std::vector<double> optimiseBasis(
				const BufferedImage<dComplex>& xy,
				const BufferedImage<double>& weight,
				const BufferedImage<double>& basis,
				const std::vector<double>& initial);

		static std::vector<double> optimiseBasis(
				const BufferedImage<dComplex>& xy,
				const BufferedImage<Tensor2x2<double>>& A,
				const BufferedImage<double>& basis,
				const std::vector<double>& initial);


		struct OldCtfBasis
		{
			double defocusU, defocusV, astigAzimuth_deg, Q0, Cs, phaseShift_deg;
		};
		
		static OldCtfBasis paramsToCtf(
				const std::vector<double> coeffs, double kV);
		
		static std::vector<double> ctfToParams(const CTF& ctf);
		
		static void testConversion();
		
		
		
		class BasisOptimisation : public Optimization
		{
			public:
		
				BasisOptimisation(
						const BufferedImage<dComplex>& xy,
						const BufferedImage<double>& weight,
						const BufferedImage<double>& basis,
						bool L1 = false);
		
				double f(const std::vector<double>& x, void* tempStorage) const;
		
				void* allocateTempStorage() const;
				void deallocateTempStorage(void* ts);
		
			private:
		
				const int w, h, cc;
				const BufferedImage<dComplex>& xy;
				const BufferedImage<double>& weight;
				const BufferedImage<double>& basis;
				const bool L1;
		};

		class AnisoBasisOptimisation : public Optimization
		{
			public:

				AnisoBasisOptimisation(
						const BufferedImage<dComplex>& xy,
						const BufferedImage<Tensor2x2<double>>& A,
						const BufferedImage<double>& basis,
						bool L1 = false);

				double f(const std::vector<double>& x, void* tempStorage) const;

				void* allocateTempStorage() const;
				void deallocateTempStorage(void* ts);

			private:

				const int w, h, cc;
				const BufferedImage<dComplex>& xy;
				const BufferedImage<Tensor2x2<double>>& A;
				const BufferedImage<double>& basis;
				const bool L1;
		};

		class MultiAnisoBasisOptimisation : public FastDifferentiableOptimization
		{
			public:

				MultiAnisoBasisOptimisation(
						const BufferedImage<dComplex>& xy,
						const BufferedImage<Tensor2x2<double>>& A,
						const BufferedImage<double>& basis,
						double lambda,
						bool L1 = false);

				double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const;
				double f(const std::vector<double>& x, void* tempStorage) const;

				void* allocateTempStorage() const;
				void deallocateTempStorage(void* ts);

			private:

				const int w, h, fc, cc;
				const double lambda;
				const BufferedImage<dComplex>& xy;
				const BufferedImage<Tensor2x2<double>>& A;
				const BufferedImage<double>& basis;
				const bool L1;
		};
};

#endif
