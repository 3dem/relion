#ifndef CUDA_HELPER_FUNCTIONS_CUH_
#define CUDA_HELPER_FUNCTIONS_CUH_
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_projector.cuh"

/*
 * This assisting function goes over the orientations determined as significant for this image, and checks
 * which translations should be included in the list of those which differences will be calculated for.
 *
 * Any contiguous translations with a shared orientation are grouped together into a "job" which is supplied
 * to the difference kernel. If there are more contiguous translations than the specified PROJDIFF_CHUNK_SIZE,
 * these are split into separate jobs, to increase paralllelism at the cost of redundant memory reads.
 */
long int makeJobsForDiff2Fine( OptimisationParamters &op,  SamplingParameters &sp,
										  long int orientation_num, long int translation_num,
					 	 	 	 	 	  ProjectionParams &FineProjectionData,
										  std::vector< long unsigned > &iover_transes,
										  std::vector< long unsigned > &ihiddens,
										  long int nr_over_orient, long int nr_over_trans, int ipart,
										  IndexedDataArray &FPW, // FPW=FinePassWeights
										  IndexedDataArrayMask &dataMask);

int  makeJobsForCollect(IndexedDataArray &FPW, IndexedDataArrayMask &dataMask);

XFLOAT thrustGetMinVal(CudaGlobalPtr<XFLOAT> &diff2s);

static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

dim3 splitCudaBlocks(long int block_num, bool doForceEven);


void mapWeights(unsigned long orientation_start,
		CudaGlobalPtr<XFLOAT> &mapped_weights,
		unsigned orientation_num, unsigned long idxArr_start, unsigned long idxArr_end,
		unsigned translation_num,
		CudaGlobalPtr <XFLOAT> &weights,
		CudaGlobalPtr <long unsigned> &rot_idx,
		CudaGlobalPtr <long unsigned> &trans_idx,
		HealpixSampling &sampling, long int ipart,
		std::vector< long unsigned > &iover_transes, std::vector< long unsigned > &ihiddens,
		std::vector< long unsigned > &iorientclasses, std::vector< long unsigned > &iover_rots,
		MultidimArray<XFLOAT> &Mweight, unsigned long current_oversampling, unsigned long nr_trans);



long unsigned imageTranslation(
		CudaGlobalPtr<XFLOAT> &Fimgs_real, CudaGlobalPtr<XFLOAT> &Fimgs_imag,
		CudaGlobalPtr<XFLOAT> &Fimgs_nomask_real, CudaGlobalPtr<XFLOAT> &Fimgs_nomask_imag,
		long int itrans_min, long int itrans_max, int adaptive_oversampling , HealpixSampling &sampling,
		std::vector<double> &oversampled_translations_x, std::vector<double> &oversampled_translations_y, std::vector<double> &oversampled_translations_z,
		unsigned long nr_oversampled_trans, std::vector<MultidimArray<Complex> > &global_fftshifts_ab_current, std::vector<MultidimArray<Complex> > &global_fftshifts_ab2_current,
		MultidimArray<Complex > &local_Fimgs_shifted, MultidimArray<Complex > &local_Fimgs_shifted_nomask,
		std::vector< long unsigned > &iover_transes, std::vector< long unsigned > &itranses, std::vector< long unsigned > &ihiddens,
		unsigned image_size);


void generateEulerMatrices(
		XFLOAT padding_factor,
		ProjectionParams ProjectionData,
		CudaGlobalPtr<XFLOAT> &eulers,
		bool inverse);

long unsigned generateProjectionSetup(
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		bool coarse,
		unsigned iclass,
		ProjectionParams &ProjectionData);


void runWavgKernel(
		CudaProjectorKernel &projector,
		CudaGlobalPtr<XFLOAT> &eulers,
		CudaGlobalPtr<XFLOAT> &Fimgs_real,
	    CudaGlobalPtr<XFLOAT> &Fimgs_imag,
	    CudaGlobalPtr<XFLOAT> &Fimgs_nomask_real,
 	    CudaGlobalPtr<XFLOAT> &Fimgs_nomask_imag,
 	    CudaGlobalPtr<XFLOAT> &sorted_weights,
 	    CudaGlobalPtr<XFLOAT> &ctfs,
 	    CudaGlobalPtr<XFLOAT> &Minvsigma2s,
 	    CudaGlobalPtr<XFLOAT> &wdiff2s_parts,
 	    CudaGlobalPtr<XFLOAT> &wavgs_real,
	    CudaGlobalPtr<XFLOAT> &wavgs_imag,
	    CudaGlobalPtr<XFLOAT> &Fweights,
	    OptimisationParamters &op,
	    MlOptimiser *baseMLO,
	    long unsigned orientation_num,
	    long unsigned translation_num,
	    unsigned image_size,
	    long int ipart,
	    int group_id,
	    int exp_iclass);

void runDiff2KernelCoarse(
		CudaProjectorKernel &projector,
		CudaGlobalPtr<XFLOAT > &gpuMinvsigma2,
		CudaGlobalPtr<XFLOAT> &Fimgs_real,
		CudaGlobalPtr<XFLOAT> &Fimgs_imag,
		CudaGlobalPtr<XFLOAT> &eulers,
		CudaGlobalPtr<XFLOAT> &diff2s,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		unsigned image_size,
		int ipart,
		int group_id,
		int exp_iclass);

void runDiff2KernelFine(
		CudaProjectorKernel &projector,
		CudaGlobalPtr<XFLOAT > &gpuMinvsigma2,
		CudaGlobalPtr<XFLOAT> &Fimgs_real,
		CudaGlobalPtr<XFLOAT> &Fimgs_imag,
		CudaGlobalPtr<XFLOAT> &eulers,
		CudaGlobalPtr<long unsigned> &rot_id,
		CudaGlobalPtr<long unsigned> &rot_idx,
		CudaGlobalPtr<long unsigned> &trans_idx,
		CudaGlobalPtr<long unsigned> &job_idx,
		CudaGlobalPtr<long unsigned> &job_num,
		CudaGlobalPtr<XFLOAT> &diff2s,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		long unsigned significant_num,
		unsigned image_size,
		int ipart,
		int group_id,
		int exp_iclass);


#endif /* CUDA_HELPER_FUNCTIONS_CUH_ */
