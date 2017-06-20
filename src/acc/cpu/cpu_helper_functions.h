#ifndef CPU_HELPER_FUNCTIONS_H_
#define CPU_HELPER_FUNCTIONS_H_

#include "src/acc/cpu/cpu_ml_optimiser.h"
#include "src/acc/cpu/cpu_projector.h"
#include "src/acc/cpu/cpu_benchmark_utils.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <signal.h>
#include "src/complex.h"
#include "src/parallel.h"
//#include <parallel_for.h>
//#include <queuing_mutex.h>

namespace CpuKernels
{
/*
 * This assisting function goes over the orientations determined as significant for this image, and checks
 * which translations should be included in the list of those which differences will be calculated for.
 *
 * Any contiguous translations with a shared orientation are grouped together into a "job" which is supplied
 * to the difference kernel. If there are more contiguous translations than the specified "chunk" number,
 * these are split into separate jobs, to increase parallelism at the cost of redundant memory reads.
 */
long int makeJobsForDiff2Fine(
		OptimisationParamters &op,  SamplingParameters &sp,
		long int orientation_num, long int translation_num,
		ProjectionParams &FineProjectionData,
		std::vector< long unsigned > &iover_transes,
		std::vector< long unsigned > &ihiddens,
		long int nr_over_orient, long int nr_over_trans, int ipart,
		IndexedDataArray &FPW, // FPW=FinePassWeights
		IndexedDataArrayMask &dataMask,
		int chunk);

/*
 * This assisting function goes over the weight-array and groups all weights with shared
 * orientations into 'jobs' which are fed into the collect-kenrel, which reduces all translations
 * with computed differences into a reduced object to be back-projected.
 */
int  makeJobsForCollect(IndexedDataArray &FPW, IndexedDataArrayMask &dataMask, unsigned long NewJobNum); // FPW=FinePassWeights

/*
 * Maps weights to a decoupled indexing of translations and orientations
 */
void mapWeights(
		unsigned long orientation_start,
		XFLOAT *mapped_weights,
		unsigned orientation_num,
		unsigned long idxArr_start,
		unsigned long idxArr_end,
		unsigned translation_num,
		XFLOAT *weights,
		long unsigned *rot_idx,
		long unsigned *trans_idx,
		unsigned long current_oversampling);

void buildCorrImage(MlOptimiser *baseMLO, OptimisationParamters &op, std::vector<XFLOAT> &corr_img, long int ipart, long int group_id);

void generateEulerMatrices(
		XFLOAT padding_factor,
		ProjectionParams &ProjectionData,
		XFLOAT *eulers,
		bool inverse);

long unsigned generateProjectionSetupFine(
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		unsigned iclass,
		ProjectionParams &ProjectionData);

void runWavgKernel(
		CpuProjectorKernel &projector,
		XFLOAT * RESTRICT eulers,
		XFLOAT * RESTRICT Fimgs_real,
		XFLOAT * RESTRICT Fimgs_imag,
		XFLOAT * RESTRICT trans_x,
		XFLOAT * RESTRICT trans_y,
		XFLOAT * RESTRICT trans_z,
		XFLOAT * RESTRICT sorted_weights,
		XFLOAT * RESTRICT ctfs,
		XFLOAT * RESTRICT wdiff2s_parts,
		XFLOAT * RESTRICT wdiff2s_AA,
		XFLOAT * RESTRICT wdiff2s_XA,
		OptimisationParamters &op,
		long unsigned orientation_num,
		long unsigned translation_num,
		unsigned image_size,
		long int ipart,
		int group_id,
		int exp_iclass,
		XFLOAT part_scale,
		bool refs_are_ctf_corrected,
		bool data_is_3D);

template< typename T>
void initValue(T *data,  size_t Size, T value)
{
	for(size_t i=0; i<Size; i++)
		data[i] = value;
}

#define WEIGHT_MAP_BLOCK_SIZE 512
void tbb_kernel_allweights_to_mweights(
		unsigned long * d_iorient,
		XFLOAT * d_allweights,
		XFLOAT * d_mweights,
		unsigned long orientation_num,
		unsigned long translation_num
		);

void mapAllWeightsToMweights(
		unsigned long * d_iorient, //projectorPlan.iorientclasses
		XFLOAT * d_allweights, //allWeights
		XFLOAT * d_mweights, //Mweight
		unsigned long orientation_num, //projectorPlan.orientation_num
		unsigned long translation_num //translation_num
		);

#define OVER_THRESHOLD_BLOCK_SIZE 512
template< typename T>
void tbb_kernel_array_over_threshold(
		int blockIdx_x, int threadIdx_x, 
		std::vector<T>    &data,
		std::vector<bool> &passed,
		T                  threshold,
		size_t             size)
{
	size_t idx = blockIdx_x * OVER_THRESHOLD_BLOCK_SIZE + threadIdx_x;
	if (idx < size)
	{
		if (data[idx] >= threshold)
			passed[idx] = true;
		else
			passed[idx] = false;
	}
}

template< typename T>
void arrayOverThreshold(std::vector<T> &data, std::vector<bool> &passed, T threshold)
{
	int grid_size = ceil((float)data.size()/(float)OVER_THRESHOLD_BLOCK_SIZE);
	for(int i=0; i<grid_size; i++) {
		for(int j=0; j<OVER_THRESHOLD_BLOCK_SIZE; j++)
			tbb_kernel_array_over_threshold<T>(
						i, j,
						data,
						passed,
						threshold,
						data.size());
	}
}

#define FIND_IN_CUMULATIVE_BLOCK_SIZE 512
template< typename T>
void tbb_kernel_find_threshold_idx_in_cumulative(
		int blockIdx_x, int threadIdx_x, 
		T *data,
		T threshold,
		size_t size_m1, //data size minus 1
		size_t *idx)
{
	size_t i = blockIdx_x * FIND_IN_CUMULATIVE_BLOCK_SIZE + threadIdx_x;
	if (i < size_m1 && data[i] <= threshold && threshold < data[i+1])
		idx[0] = i+1;
}

size_t findThresholdIdxInCumulativeSum(std::vector<XFLOAT> &data, XFLOAT threshold);

void runDiff2KernelCoarse(
		CpuProjectorKernel &projector,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
		XFLOAT *corr_img,
		XFLOAT *Fimg_real,
		XFLOAT *Fimg_imag,
		XFLOAT *d_eulers,
		XFLOAT *diff2s,
		XFLOAT local_sqrtXi2,
		long unsigned orientation_num,
		int translation_num,
		int image_size,
		bool do_CC,
		bool data_is_3D);

void runDiff2KernelFine(
		CpuProjectorKernel &projector,
		XFLOAT *corr_img,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
		XFLOAT *eulers,
		long unsigned *rot_id,
		long unsigned *rot_idx,
		long unsigned *trans_idx,
		long unsigned *job_idx,
		long unsigned *job_num,
		XFLOAT *diff2s,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned translation_num,
		unsigned image_size,
		int ipart,
		int exp_iclass,
		long unsigned job_num_count,
		bool do_CC,
		bool data_is_3D);

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
void tbb_kernel_window_fourier_transform(
					int blockIdx_x,
					int blockIdx_y,  
					int threadIdx_x, 
		XFLOAT *g_in_real,
		XFLOAT *g_in_imag,
		XFLOAT *g_out_real,
		XFLOAT *g_out_imag,
		unsigned iX, unsigned iY, unsigned iZ, unsigned iYX, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ, unsigned oYX, //Output dimensions
		unsigned max_idx,
		unsigned max_r2 = 0
		)
{
	unsigned n = threadIdx_x + WINDOW_FT_BLOCK_SIZE * blockIdx_x;
	long int image_offset = oX*oY*oZ*blockIdx_y;
	if (n >= max_idx) return;

	int k, i, kp, ip, jp;

	if (check_max_r2)
	{
		k = n / (iX * iY);
		i = (n % (iX * iY)) / iX;

		kp = k < iX ? k : k - iZ;
		ip = i < iX ? i : i - iY;
		jp = n % iX;

		if (kp*kp + ip*ip + jp*jp > max_r2)
			return;
	}
	else
	{
		k = n / (oX * oY);
		i = (n % (oX * oY)) / oX;

		kp = k < oX ? k : k - oZ;
		ip = i < oX ? i : i - oY;
		jp = n % oX;
	}

	g_out_real[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp + image_offset] = g_in_real[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp + image_offset];
	g_out_imag[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp + image_offset] = g_in_imag[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp + image_offset];
}

void runCollect2jobs(	int      grid_dim,
						XFLOAT * oo_otrans_x,          // otrans-size -> make const
						XFLOAT * oo_otrans_y,          // otrans-size -> make const
						XFLOAT * oo_otrans_z,          // otrans-size -> make const
						XFLOAT * myp_oo_otrans_x2y2z2, // otrans-size -> make const
						XFLOAT * weights,
						XFLOAT significant_weight,    // TODO Put in const
						XFLOAT sum_weight,    		  // TODO Put in const
						unsigned long nr_trans,
						unsigned long oversampled_trans,
						unsigned long oversampled_rot,
						int oversamples,
						bool skip_rots,
						XFLOAT * p_weights,
						XFLOAT * p_thr_wsum_prior_offsetx_class,
						XFLOAT * p_thr_wsum_prior_offsety_class,
						XFLOAT * p_thr_wsum_prior_offsetz_class,
						XFLOAT * p_thr_wsum_sigma2_offset,
						size_t * rot_idx,
						size_t * trans_idx,
						size_t * jobOrigin,
						size_t * jobExtent,
						bool data_is_3D
						);

void windowFourierTransform2(
		XFLOAT *d_in_real,
		XFLOAT *d_in_imag,
		XFLOAT *d_out_real,
		XFLOAT *d_out_imag,
		unsigned iX, unsigned iY, unsigned iZ, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ  //Output dimensions
					);

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
void tbb_kernel_window_fourier_transform(
		int blockIdx_x,
		int blockIdx_y,  
		int threadIdx_x,
		CUDACOMPLEX *g_in,
		CUDACOMPLEX *g_out,
		size_t iX, size_t iY, size_t iZ, size_t iYX, //Input dimensions
		size_t oX, size_t oY, size_t oZ, size_t oYX, //Output dimensions
		size_t max_idx,
		size_t max_r2 = 0
		)
{
	size_t n = threadIdx_x + WINDOW_FT_BLOCK_SIZE * blockIdx_x;
	size_t oOFF = oX*oY*oZ*blockIdx_y;
	size_t iOFF = iX*iY*iZ*blockIdx_y;
	if (n >= max_idx) return;

	long int k, i, kp, ip, jp;

	if (check_max_r2)
	{
		k = n / (iX * iY);
		i = (n % (iX * iY)) / iX;

		kp = k < iX ? k : k - iZ;
		ip = i < iX ? i : i - iY;
		jp = n % iX;

		if (kp*kp + ip*ip + jp*jp > max_r2)
			return;
	}
	else
	{
		k = n / (oX * oY);
		i = (n % (oX * oY)) / oX;

		kp = k < oX ? k : k - oZ;
		ip = i < oX ? i : i - oY;
		jp = n % oX;
	}

	long int  in_idx = (kp < 0 ? kp + iZ : kp) * iYX + (ip < 0 ? ip + iY : ip)*iX + jp;
	long int out_idx = (kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp;
	g_out[out_idx + oOFF] =  g_in[in_idx + iOFF];
}

void windowFourierTransform2(
		std::vector<CUDACOMPLEX > &d_in,
		std::vector<CUDACOMPLEX > &d_out,
		size_t iX, size_t iY, size_t iZ, //Input dimensions
		size_t oX, size_t oY, size_t oZ,  //Output dimensions
		size_t Npsi = 1,
		size_t pos = 0);


void selfApplyBeamTilt2(MultidimArray<Complex > &Fimg, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
		RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size);

template <typename T>
void runCenterFFT(MultidimArray< T >& v, bool forward)
{
	std::vector<XFLOAT >  img_in (v.nzyxdim);   // with original data pointer
//	std::vector<XFLOAT >  img_aux(v.nzyxdim, allocator);   // temporary holder

	for (unsigned i = 0; i < v.nzyxdim; i ++)
		img_in[i] = (XFLOAT) v.data[i];

	if ( v.getDim() == 1 )
	{
		std::cerr << "CenterFFT on gpu reverts to cpu for dim!=2 (now dim=1)" <<std::endl;
		// 1D
		MultidimArray< T > aux;
		int l, shift;

		l = XSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		// Shift the input in an auxiliar vector
		for (int i = 0; i < l; i++)
		{
			int ip = i + shift;

			if (ip < 0)
				ip += l;
			else if (ip >= l)
				ip -= l;

			aux(ip) = DIRECT_A1D_ELEM(v, i);
		}

		// Copy the vector
		for (int i = 0; i < l; i++)
			DIRECT_A1D_ELEM(v, i) = DIRECT_A1D_ELEM(aux, i);
	}
	else if ( v.getDim() == 2 )
	{
		// 2D
		//std::cerr << "CenterFFT on gpu with dim=2!" <<std::endl;

		long int xshift = (int)(XSIZE(v) / 2);
		long int yshift = (int)(YSIZE(v) / 2);

		if (!forward)
		{
			xshift = -xshift;
			yshift = -yshift;
		}

		//TODO-SFU - batch size assumed 1, so second index is forced to 0
		// Is this right???
		int  dim = ceilf((float)(v.nzyxdim/(float)(2*CFTT_BLOCK_SIZE)));
		for(int i=0; i<dim; i++) {
			for(int k=0; k<CFTT_BLOCK_SIZE; k++)
				tbb_kernel_centerFFT_2D(i, 0, k, img_in,
				v.nzyxdim, XSIZE(v), YSIZE(v), xshift, yshift);
		}

		for (unsigned i = 0; i < v.nzyxdim; i ++)
			v.data[i] = (T) img_in[i];

	}
	else if ( v.getDim() == 3 )
	{
		std::cerr << "CenterFFT on gpu reverts to cpu for dim!=2 (now dim=3)" <<std::endl;
		// 3D
		MultidimArray< T > aux;
		int l, shift;

		// Shift in the X direction
		l = XSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		for (int k = 0; k < ZSIZE(v); k++)
			for (int i = 0; i < YSIZE(v); i++)
			{
				// Shift the input in an auxiliar vector
				for (int j = 0; j < l; j++)
				{
					int jp = j + shift;

					if (jp < 0)
						jp += l;
					else if (jp >= l)
						jp -= l;

					aux(jp) = DIRECT_A3D_ELEM(v, k, i, j);
				}

				// Copy the vector
				for (int j = 0; j < l; j++)
					DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, j);
			}

		// Shift in the Y direction
		l = YSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		for (int k = 0; k < ZSIZE(v); k++)
			for (int j = 0; j < XSIZE(v); j++)
			{
				// Shift the input in an auxiliar vector
				for (int i = 0; i < l; i++)
				{
					int ip = i + shift;

					if (ip < 0)
						ip += l;
					else if (ip >= l)
						ip -= l;

					aux(ip) = DIRECT_A3D_ELEM(v, k, i, j);
				}

				// Copy the vector
				for (int i = 0; i < l; i++)
					DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, i);
			}

		// Shift in the Z direction
		l = ZSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		for (int i = 0; i < YSIZE(v); i++)
			for (int j = 0; j < XSIZE(v); j++)
			{
				// Shift the input in an auxiliar vector
				for (int k = 0; k < l; k++)
				{
					int kp = k + shift;
					if (kp < 0)
						kp += l;
					else if (kp >= l)
						kp -= l;

					aux(kp) = DIRECT_A3D_ELEM(v, k, i, j);
				}

				// Copy the vector
				for (int k = 0; k < l; k++)
					DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, k);
			}
	}
	else
	{
		v.printShape();
		REPORT_ERROR("CenterFFT ERROR: Dimension should be 1, 2 or 3");
	}
}


template <typename T>
void runCenterFFT( std::vector< T > &img_in,
			 int  xSize,
			 int  ySize,
			 bool forward,
			 int  batchSize = 1)
{
	int xshift = (xSize / 2);
	int yshift = (ySize / 2);

	if (!forward)
	{
		xshift = -xshift;
		yshift = -yshift;
	}

	int blocks = ceilf((float)((xSize*ySize)/(float)(2*CFTT_BLOCK_SIZE)));
	for(int i=0; i<blocks; i++) {
		for(int j=0; j<batchSize; j++)
			for(int k=0; k<CFTT_BLOCK_SIZE; k++)
				tbb_kernel_centerFFT_2D(i, j, k,
					&img_in[0],
					 xSize*ySize,
					 xSize,
					 ySize,
					 xshift,
					 yshift);
	};

//	HANDLE_ERROR(cudaStreamSynchronize(0));
//	img_aux.cp_on_device(img_in.d_ptr); //update input image with centered kernel-output.


}

template <typename T>
void runCenterFFT(std::vector< T > &img_in,
			int xSize,
			int ySize,
			int zSize,
			bool forward,
			int  batchSize = 1)
{
//	CudaGlobalPtr<XFLOAT >  img_aux(img_in.h_ptr, img_in.size, allocator);   // temporary holder
//	img_aux.device_alloc();

	if(zSize>1)
	{
		int xshift = (xSize / 2);
		int yshift = (ySize / 2);
		int zshift = (ySize / 2);

		if (!forward)
		{
			xshift = -xshift;
			yshift = -yshift;
			zshift = -zshift;
		}

		int block =ceilf((float)((xSize*ySize*zSize)/(float)(2*CFTT_BLOCK_SIZE)));
		for(int i=0; i<block; i++){
			for(int j=0; j<batchSize; j++)
				for(int k=0; k<CFTT_BLOCK_SIZE; k++)
					tbb_kernel_centerFFT_3D(
						i,j,k,
						&img_in[0],
						xSize*ySize*zSize,
						xSize,
						ySize,
						zSize,
						xshift,
						yshift,
						zshift);
		}

	}
	else
	{
		int xshift = (xSize / 2);
		int yshift = (ySize / 2);

		if (!forward)
		{
			xshift = -xshift;
			yshift = -yshift;
		}

		int block = ceilf((float)((xSize*ySize)/(float)(2*CFTT_BLOCK_SIZE)));
		for(int i=0; i<block; i++){
			for(int j=0; j<batchSize; j++)
				for(int k=0; k<CFTT_BLOCK_SIZE; k++)
					tbb_kernel_centerFFT_2D(
						i, j, k,
						&img_in[0],
						xSize*ySize,
						xSize,
						ySize,
						xshift,
						yshift);
		}
	}

}

template <typename T>
void lowPassFilterMapGPU(
		std::vector< T > &img_in,
		size_t Zdim,
		size_t Ydim,
		size_t Xdim,
		long int ori_size,
		RFLOAT lowpass,
		RFLOAT highpass,
		RFLOAT angpix,
		int filter_edge_width,
		bool do_highpass)
{
	// High or low?
	RFLOAT passLimit = (do_highpass ? highpass : lowpass);

	// Which resolution shell is the filter?
	int ires_filter = ROUND((ori_size * angpix)/passLimit);
	int filter_edge_halfwidth = filter_edge_width / 2;

	// Soft-edge: from 1 shell less to one shell more:
	XFLOAT edge_low = XMIPP_MAX(0., (ires_filter - filter_edge_halfwidth) / (RFLOAT)ori_size); // in 1/pix
	XFLOAT edge_high = XMIPP_MIN(Xdim, (ires_filter + filter_edge_halfwidth) / (RFLOAT)ori_size); // in 1/pix
	XFLOAT edge_width = edge_high - edge_low;

	int blocks =ceilf( (float)(Xdim*Ydim*Zdim)/ (float)(CFTT_BLOCK_SIZE) );
	if (do_highpass)
	{
	   for(int i=0; i<blocks; i++) {
			for(int j=0; j<CFTT_BLOCK_SIZE; j++)
			   tbb_kernel_frequencyPass<true>(i, j,
				img_in,
				ori_size,
				Xdim,
				Ydim,
				Zdim,
				edge_low,
				edge_width,
				edge_high,
				(XFLOAT)angpix,
				Xdim*Ydim*Zdim);
		}
	}
	else
	{
		for(int i=0; i<blocks; i++) {
			for(int j=0; j<CFTT_BLOCK_SIZE; j++)
				tbb_kernel_frequencyPass<false>(i, j,
						img_in,
						ori_size,
						Xdim,
						Ydim,
						Zdim,
						edge_low,
						edge_width,
						edge_high,
						(XFLOAT)angpix,
						Xdim*Ydim*Zdim);
		}
	}
}

// may need to change to parallel reduce if it becomes the bottle neck.
template <typename T>
static T getMin(T *data, size_t size)
{
	T min = data[0];
	for(int i=1; i<size; i++)
		min = data[i] < min ? data[i] : min;
	 
	return min;
}

template <typename T>
static T getSum(T *data, size_t size)
{
	T sum = data[0];
	for(int i=1; i<size; i++)
		sum += data[i];
	 
	return sum;
}

template <typename T>
static std::pair<int, T> getArgMin(T *data, size_t size)
{
	std::pair<int, T> pair;
	pair.first = 0;
	pair.second = data[0];
	
	for(int i=1; i<size; i++)
		if( data[i] < pair.second) {
			pair.first = i;
			pair.second = data[i];
		}
		
	return pair;
}

template <typename T>
static std::pair<int, T> getArgMax(T *data, size_t size)
{
	std::pair<int, T> pair;
	pair.first = 0;
	pair.second = data[0];
	
	for(int i=1; i<size; i++)
		if( data[i] > pair.second) {
			pair.first = i;
			pair.second = data[i];
		}
		
	return pair;
}
} // Namespace CpuKernels

#endif //CPU_HELPER_FUNCTIONS_H_

