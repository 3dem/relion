/*
 * cuda_autopicker.h
 *
 *  Created on: Dec 3, 2015
 *      Author: bjornf
 */

#ifndef CUDA_AUTOPICKER_H_
#define CUDA_AUTOPICKER_H_

#include "src/mpi.h"
#include "src/autopicker.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/projector.h"
#include "src/gpu_utils/cuda_projector.h"

#include <stack>

#ifdef CUDA_DOUBLE_PRECISION
#define XFLOAT double
#else
#define XFLOAT float
#endif

class AutoPickerCuda
{
public:

	AutoPicker *basePckr;

	CudaCustomAllocator *allocator;

	std::vector< CudaProjector > cudaProjectors;

   //Class streams ( for concurrent scheduling of class-specific kernels)
	std::vector< cudaStream_t > classStreams;

	int device_id;

	//MlDeviceBundle *devBundle;

	AutoPickerCuda(AutoPicker *basePicker, int dev_id);

	void setupProjectors();

	void run();

	void autoPickOneMicrograph(FileName &fn_mic);

	~AutoPickerCuda()
	{
		for (int i = 0; i < classStreams.size(); i++)
			HANDLE_ERROR(cudaStreamDestroy(classStreams[i]));
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

#endif /* CUDA_AUTOPICKER_H_ */
