#ifndef CUDA_SKUNKS_CUH_
#define CUDA_SKUNKS_CUH_

#include "src/projector.h"
#include "src/multidim_array.h"
#include "src/fftw.h"

void computeFourierTransformMap(Projector *P, MultidimArray<RFLOAT> &vol_in, MultidimArray<RFLOAT> &power_spectrum, int current_size = -1, int nr_threads = 1, bool do_gridding = true)
{
	MultidimArray<RFLOAT> Mpad;
	MultidimArray<Complex > Faux;
	FourierTransformer transformer;
	RFLOAT normfft;

	// Size of padded real-space volume
	int padoridim = P->padding_factor * P->ori_size;

	// Initialize data array of the oversampled transform
	P->ref_dim = vol_in.getDim();

	// Make Mpad
	switch (P->ref_dim)
	{
	case 2:
	   Mpad.initZeros(padoridim, padoridim);
	   normfft = (RFLOAT)(P->padding_factor * P->padding_factor);
	   break;
	case 3:
	   Mpad.initZeros(padoridim, padoridim, padoridim);
	   if (P->data_dim ==3)
		   normfft = (RFLOAT)(P->padding_factor * P->padding_factor * P->padding_factor);
	   else
		   normfft = (RFLOAT)(P->padding_factor * P->padding_factor * P->padding_factor * P->ori_size);
	   break;
	default:
	   REPORT_ERROR("Projector::computeFourierTransformMap%%ERROR: Dimension of the data array should be 2 or 3");
	}

	// First do a gridding pre-correction on the real-space map:
	// Divide by the inverse Fourier transform of the interpolator in Fourier-space
	// 10feb11: at least in 2D case, this seems to be the wrong thing to do!!!
	// TODO: check what is best for subtomo!
	if (do_gridding)// && data_dim != 3)
		P->griddingCorrect(vol_in);

	// Pad translated map with zeros
	vol_in.setXmippOrigin();
	Mpad.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in) // This will also work for 2D
		A3D_ELEM(Mpad, k, i, j) = A3D_ELEM(vol_in, k, i, j);

	// Translate padded map to put origin of FT in the center
	CenterFFT(Mpad, true);

	// Calculate the oversampled Fourier transform
	transformer.FourierTransform(Mpad, Faux, false);

	// Free memory: Mpad no longer needed
	Mpad.clear();

	// Resize data array to the right size and initialise to zero
	P->initZeros(current_size);

	// Fill data only for those points with distance to origin less than max_r
	// (other points will be zero because of initZeros() call above
	// Also calculate radial power spectrum
	power_spectrum.initZeros(P->ori_size / 2 + 1);
	MultidimArray<RFLOAT> counter(power_spectrum);
	counter.initZeros();

	int max_r2 = P->r_max * P->r_max * P->padding_factor * P->padding_factor;
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) // This will also work for 2D
	{
		int r2 = kp*kp + ip*ip + jp*jp;
		// The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
		if (r2 <= max_r2)
		{
			// Set data array
			A3D_ELEM(P->data, kp, ip, jp) = DIRECT_A3D_ELEM(Faux, k, i, j) * normfft;

			// Calculate power spectrum
			int ires = ROUND( sqrt((RFLOAT)r2) / P->padding_factor );
			// Factor two because of two-dimensionality of the complex plane
			DIRECT_A1D_ELEM(power_spectrum, ires) += norm(A3D_ELEM(P->data, kp, ip, jp)) / 2.;
			DIRECT_A1D_ELEM(counter, ires) += 1.;
		}
	}

	// Calculate radial average of power spectrum
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(power_spectrum)
	{
		if (DIRECT_A1D_ELEM(counter, i) < 1.)
			DIRECT_A1D_ELEM(power_spectrum, i) = 0.;
		else
			DIRECT_A1D_ELEM(power_spectrum, i) /= DIRECT_A1D_ELEM(counter, i);
	}
}

#endif //CUDA_SKUNKS_CUH_

