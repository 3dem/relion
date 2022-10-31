/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
#ifndef HIP_AUTOPICKER_H_
#define HIP_AUTOPICKER_H_

#include "src/mpi.h"
#include "src/autopicker.h"
#include "src/autopicker_mpi.h"
#include "src/projector.h"
#include "src/complex.h"
#include "src/image.h"

#include "src/acc/hip/hip_mem_utils.h"
#include "src/acc/acc_projector.h"
#include "src/acc/hip/hip_settings.h"
#include "src/acc/hip/hip_fft.h"
#include "src/acc/hip/hip_benchmark_utils.h"

#include <stack>

#ifdef ACC_DOUBLE_PRECISION
#define XFLOAT double
#else
#define XFLOAT float
#endif

class AutoPickerHip
{
private:

	MpiNode *node;

public:

	AutoPicker *basePckr;

	HipCustomAllocator *allocator;
	HipFFT micTransformer;
	HipFFT hipTransformer1;
	HipFFT hipTransformer2;

	std::vector< AccProjector > projectors;

   //Class streams ( for concurrent scheduling of class-specific kernels)
	std::vector< hipStream_t > classStreams;

	int device_id;

	bool have_warned_batching;

	//MlDeviceBundle *devBundle;

#ifdef TIMING_FILES
	relion_timer timer;
#endif

	AutoPickerHip(AutoPicker    *basePicker, const char * timing_fnm);
	AutoPickerHip(AutoPickerMpi *basePicker, const char * timing_fnm);

	void setupProjectors();

	void run();

	void autoPickOneMicrograph(FileName &fn_mic, long int imic);

	void calculateStddevAndMeanUnderMask(AccPtr< ACCCOMPLEX > &d_Fmic,
			AccPtr< ACCCOMPLEX > &d_Fmic2,
			AccPtr< ACCCOMPLEX > &d_Fmsk,
			int nr_nonzero_pixels_mask, AccPtr< XFLOAT > &d_Mstddev,
			AccPtr< XFLOAT > &d_Mmean,
			size_t x, size_t y, size_t mic_size, size_t workSize);

	~AutoPickerHip()
	{
		for (int i = 0; i < classStreams.size(); i++)
			HANDLE_ERROR(hipStreamDestroy(classStreams[i]));
	}

//private:

//	// Uses Roseman2003 formulae to calculate stddev under the mask through FFTs
//	// The FFTs of the micrograph (Fmic), micrograph-squared (Fmic2) and the mask (Fmsk) need to be provided at downsize_mic
//	// The putput (Mstddev) will be at (binned) micrograph_size
//	void calculateStddevAndMeanUnderMask(const MultidimArray<Complex > &Fmic, const MultidimArray<Complex > &Fmic2,
//			MultidimArray<Complex > &Fmsk, int nr_nonzero_pixels_mask, MultidimArray<RFLOAT> &Mstddev, MultidimArray<RFLOAT> &Mmean);
//
//	// Peak search for all pixels above a given threshold in the map
//	void peakSearch(const MultidimArray<RFLOAT> &Mccf, const MultidimArray<RFLOAT> &Mpsi, const MultidimArray<RFLOAT> &Mstddev, int iref, int skip_side, std::vector<Peak> &peaks);
//
//	// Now prune the coordinates: within min_particle_distance: all peaks are the same cluster
//	// From each cluster, take the single peaks with the highest ccf
//	// If then, there is another peaks at a distance of at least min_particle_distance: take that one as well, and so forth...
//	void prunePeakClusters(std::vector<Peak> &peaks, int min_distance);
//
//
//	// Only keep those peaks that are at the given distance apart from each other
//	void removeTooCloselyNeighbouringPeaks(std::vector<Peak> &peaks, int min_distance);

};

#endif /* HIP_AUTOPICKER_H_ */
