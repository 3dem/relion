#ifndef HELPER_GPU_KERNELS_H_
#define HELPER_GPU_KERNELS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <cassert>
#include <sycl/sycl.hpp>
#ifdef USE_ONEDPL
// Needs oneDPL (https://github.com/oneapi-src/oneDPL)
 #include <oneapi/dpl/algorithm>
 #include <oneapi/dpl/execution>
 #include <oneapi/dpl/iterator>
 #include <oneapi/dpl/numeric>
 #include <oneapi/dpl/ranges>
#endif

#include "src/macros.h"
#include "src/acc/sycl/sycl_settings.h"
#include "src/acc/sycl/sycl_dev.h"
#include "src/acc/sycl/sycl_kernels/sycl_utils.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projectorkernel_impl.h"

namespace syclGpuKernels
{
using atomic_ref_dev = sycl::atomic_ref<XFLOAT, sycl::memory_order::relaxed, sycl::memory_scope::device, sycl::access::address_space::global_space>;
using atomic_ref_wg = sycl::atomic_ref<XFLOAT, sycl::memory_order::relaxed, sycl::memory_scope::work_group, sycl::access::address_space::global_space>;
using atomic_ref_sg = sycl::atomic_ref<XFLOAT, sycl::memory_order::relaxed, sycl::memory_scope::sub_group, sycl::access::address_space::global_space>;

template <typename I>
constexpr void
__group_barrier(I __item)
{
#if 0
// group_barrier is too slow!
	sycl::group_barrier(__item.get_group(), sycl::memory_scope::work_group);
#else
	__item.barrier(sycl::access::fence_space::local_space);
#endif
}

template< typename T>
void InitValue(T *data, const T& value, const size_t count, virtualSYCL *devAcc)
{
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	Q->fill<T>(data, value, count).wait_and_throw();
}

template <typename T>
static T getSumOnDevice(T *data, const size_t count, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	T ret;
 #if defined(USE_ONEDPL) && !defined(USE_LESS_ONEDPL)
	ret = dpl::reduce(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), data, data+count, static_cast<T>(0));
 #else
	{
		using namespace sycl;
		auto buf = buffer{&ret, range{1}};
		Q->submit([&](handler &cgh)
		{
			auto reductionSum = reduction(buf, cgh, plus<T>{}, property::reduction::initialize_to_identity{});
			cgh.parallel_for(range<1>{count}, reductionSum, [=](id<1> idx, auto& sum)
			{
				sum += data[idx];
			});
		}).wait_and_throw();
	}
 #endif
	return ret;
}

template <typename T>
static T getMinOnDevice(T *data, const size_t count, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	T val;
 #if defined(USE_ONEDPL) && !defined(USE_LESS_ONEDPL)
	auto iter = dpl::min_element(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), data, data+count);
	auto dist = dpl::distance(data, iter);

	Q->memcpy(&val, data+dist, sizeof(T)).wait_and_throw();
 #else
	val = std::numeric_limits<T>::max();
	{
		using namespace sycl;
		auto buf = buffer{&val, range{1}};
		Q->submit([&](handler &cgh)
		{
			auto reductionMin = reduction(buf, cgh, minimum<T>{});
			cgh.parallel_for(range<1>{count}, reductionMin, [=](id<1> idx, auto& vmin)
			{
				vmin.combine(data[idx]);
			});
		}).wait_and_throw();
	}
 #endif
	return val;
}

template <typename T>
static T getMaxOnDevice(T *data, const size_t count, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	T val;
 #if defined(USE_ONEDPL) && !defined(USE_LESS_ONEDPL)
	auto iter = dpl::max_element(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), data, data+count);
	auto dist = dpl::distance(data, iter);

	Q->memcpy(&val, data+dist, sizeof(T)).wait_and_throw();
 #else
	val = std::numeric_limits<T>::lowest();
	{
		using namespace sycl;
		auto buf = buffer{&val, range{1}};
		Q->submit([&](handler &cgh)
		{
			auto reductionMax = reduction(buf, cgh, maximum<T>{});
			cgh.parallel_for(range<1>{count}, reductionMax, [=](id<1> idx, auto& vmax)
			{
				vmax.combine(data[idx]);
			});
		}).wait_and_throw();
	}
 #endif
	return val;
}

template <typename T, typename I>
using minloc = sycl::minimum<std::pair<T, I>>;

template <typename T>
static std::pair<size_t, T> getArgMinOnDevice(T *data, const size_t count, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	std::pair<size_t, T> pair;
 #if defined(USE_ONEDPL) && !defined(USE_LESS_ONEDPL)
	auto iter = dpl::min_element(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), data, data+count);
	auto dist = dpl::distance(data, iter);

	pair.first = dist;
	Q->memcpy(&pair.second, data+dist, sizeof(T)).wait_and_throw();
 #else
	const std::pair<T,size_t> identity { std::numeric_limits<T>::max(), std::numeric_limits<size_t>::max() };
	std::pair<T,size_t> res { identity };
	{
		using namespace sycl;
		auto buf = buffer{&res, range{1}};
		Q->submit([&](handler &cgh)
		{
			auto reductionMin = reduction(buf, cgh, identity, minloc<T,size_t>{});
			cgh.parallel_for(range<1>{count}, reductionMin, [=](id<1> idx, auto& vmin)
			{
				std::pair<T,size_t> temp { data[idx], idx };
				vmin.combine(temp);
			});
		}).wait_and_throw();
	}

	pair.first = res.second;
	pair.second = res.first;
 #endif
	return pair;
}

template <typename T, typename I>
using maxloc = sycl::maximum<std::pair<T, I>>;

template <typename T>
static std::pair<size_t, T> getArgMaxOnDevice(T *data, const size_t count, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	std::pair<size_t, T> pair;
 #if defined(USE_ONEDPL) && !defined(USE_LESS_ONEDPL)
	auto iter = dpl::max_element(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), data, data+count);
	auto dist = dpl::distance(data, iter);

	pair.first = dist;
	Q->memcpy(&pair.second, data+dist, sizeof(T)).wait_and_throw();
 #else
	const std::pair<T,size_t> identity { std::numeric_limits<T>::lowest(), std::numeric_limits<size_t>::lowest() };
	std::pair<T,size_t> res { identity };
	{
		using namespace sycl;
		auto buf = buffer{&res, range{1}};
		Q->submit([&](handler &cgh)
		{
			auto reductionMax = reduction(buf, cgh, identity, maxloc<T,size_t>{});
			cgh.parallel_for(range<1>{count}, reductionMax, [=](id<1> idx, auto& vmax)
			{
				std::pair<T,size_t> temp { data[idx], idx };
				vmax.combine(temp);
			});
		}).wait_and_throw();
	}

	pair.first = res.second;
	pair.second = res.first;
 #endif
	return pair;
}

template <typename T>
static size_t filterGreaterZeroOnDevice(T *in, const size_t count, T *out, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	size_t dist;
 #if defined(USE_ONEDPL) && !defined(USE_LESS_ONEDPL)
	auto iter = dpl::copy_if(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), in, in+count, out, [=](T &v){return (v>static_cast<T>(0));});
	dist = dpl::distance(out, iter);
 #else
	dist = 0;
	{
		using namespace sycl;
		using atomic_ref_count = atomic_ref<size_t, memory_order::relaxed, memory_scope::device, access::address_space::global_space>;

		auto buf = buffer{&dist, range{1}};

		Q->submit([&](handler &cgh)
		{
			auto ncopy = accessor{buf, cgh, read_write};
			cgh.parallel_for(range<1>(count), [=](id<1> idx)
			{
				if (in[idx] > (T)0)
				{
					auto i = atomic_ref_count(ncopy[0]).fetch_add(1);
					out[i] = in[idx];
				}
			});
		}).wait_and_throw();
	}
 #endif
	return dist;
}

template <typename T>
static void sortOnDevice(T *in, const size_t count, T *out, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
 #ifdef USE_ONEDPL
	Q->memcpy(out, in, count *sizeof(T)).wait_and_throw();
	dpl::sort(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), out, out+count);
 #endif
}

template <typename T>
static void scanOnDevice(T *in, const size_t count, T *out, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
 #ifdef USE_MORE_ONEDPL
// This has conflict with DisableIndirectAccess=1
// DisableIndirectAccess=1 will cause incorrecxt result
	dpl::inclusive_scan(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), in, in+count, out);
 #else
	{
// This requires SYCL2020 for joint_inclusive_scan
		using namespace sycl;

		const std::array<int,3> maxItem = dynamic_cast<devSYCL*>(devAcc)->maxItem;
		const size_t runWorkItem = (count > maxItem[2]) ? maxItem[2] : count;
		const size_t G = (count%runWorkItem == 0) ? count/runWorkItem : count/runWorkItem + 1;
		const size_t WG = (G <= maxItem[1]) ? G : maxItem[1];
		const size_t ND = (G <= maxItem[1]) ? 1 : ((G%WG == 0) ? G/WG : G/WG + 1);

		Q->submit([&](handler& cgh)
		{
			cgh.parallel_for(nd_range<3>(range<3>(ND,WG,runWorkItem), range<3>(1,1,runWorkItem)), [=](nd_item<3> it)
			{
				auto g = it.get_group();
				joint_inclusive_scan(g, in, in+count, out, sycl::plus<T>{});
			});
		}).wait_and_throw();
	}
 #endif
}

static void sycl_kernel_exponentiate_weights_fine(
		XFLOAT *g_pdf_orientation,
		bool *g_pdf_orientation_zeros,
		XFLOAT *g_pdf_offset,
		bool *g_pdf_offset_zeros,
		XFLOAT *g_weights,
		const XFLOAT min_diff2,
		const unsigned long oversamples_orient,
		const unsigned long oversamples_trans,
		unsigned long *d_rot_id,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		const long job_num,
		virtualSYCL *devAcc)
{
	assert(job_num <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		using namespace sycl;

		cgh.parallel_for(range<1>(job_num), [=](id<1> jobid)
		{
			unsigned long pos = d_job_idx[jobid];
			// index of comparison
			unsigned long ix = d_rot_id   [pos];   // each thread gets its own orient...
			unsigned long iy = d_trans_idx[pos];   // ...and it's starting trans...
			int in = static_cast<int>(d_job_num[jobid]); // ...AND the number of translations to go through

			for (int itrans=0; itrans < in; itrans++, iy++)
			{
				unsigned long c_itrans = ( iy - (iy % oversamples_trans))/ oversamples_trans;

				if( g_weights[pos+itrans] < min_diff2 || g_pdf_orientation_zeros[ix] || g_pdf_offset_zeros[c_itrans])
					g_weights[pos+itrans] = -99e99; //large negative number
				else
					g_weights[pos+itrans] = g_pdf_orientation[ix] + g_pdf_offset[c_itrans] + min_diff2 - g_weights[pos+itrans];
			}
		});
	}).wait_and_throw();
}

 #ifdef ACC_DOUBLE_PRECISION
 # define KERNEL_EXP_VALUE (-700.0)
 #else
 # define KERNEL_EXP_VALUE (-88.0f)
 #endif
template<typename T>
static void sycl_kernel_exponentiate(
		T *g_array,
		const T add,
		const size_t count,
		virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
 #ifdef USE_MORE_ONEDPL
  #ifdef ACC_DOUBLE_PRECISION
	dpl::transform(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), g_array, g_array+count, g_array, [=](T v) { v+=add; v = (v<KERNEL_EXP_VALUE) ? static_cast<T>(0) : sycl::exp(v); return v;});
  #else
	dpl::transform(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), g_array, g_array+count, g_array, [=](T v) { v+=add; v = (v<KERNEL_EXP_VALUE) ? static_cast<T>(0) : sycl::native::exp(v); return v;});
  #endif
 #else
	Q->submit([&](sycl::handler &cgh)
	{
		using namespace sycl;

		cgh.parallel_for(range<1>(count), [=](id<1> idx)
		{
			if (idx < count)
			{
				T a = g_array[idx] + add;
				if (a < KERNEL_EXP_VALUE)
					g_array[idx] = static_cast<T>(0);
				else
  #ifdef ACC_DOUBLE_PRECISION
					g_array[idx] = sycl::exp(a);
  #else
					g_array[idx] = sycl::native::exp(a);
  #endif
			}
		});
	}).wait_and_throw();
 #endif
}

template<typename T>
static void sycl_kernel_weights_exponent_coarse(
		T *g_pdf_orientation,
		bool *g_pdf_orientation_zeros,
		T *g_pdf_offset,
		bool *g_pdf_offset_zeros,
		T *g_weights,
		const T min_diff2,
		const int nr_coarse_orient,
		const int nr_coarse_trans,
		const unsigned long max_idx,
		virtualSYCL *devAcc)
{
	assert(max_idx <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		using namespace sycl;

		cgh.parallel_for(range<1>(max_idx), [=](id<1> idx)
		{
			const int itrans = idx % nr_coarse_trans;
			const int iorient = (idx - itrans) / nr_coarse_trans;

			T diff2 = g_weights[idx];
			if( diff2 < min_diff2 || g_pdf_orientation_zeros[iorient] || g_pdf_offset_zeros[itrans])
				g_weights[idx] = -99e99; //large negative number
			else
				g_weights[idx] = g_pdf_orientation[iorient] + g_pdf_offset[itrans] + min_diff2 - diff2;
		});
	}).wait_and_throw();
}

template<bool DATA3D>
inline void collect2jobs(  int     grid_size,
					int		block_size,
					XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
					XFLOAT *g_oo_otrans_y,          // otrans-size -> make const
					XFLOAT *g_oo_otrans_z,          // otrans-size -> make const
					XFLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
					XFLOAT *g_i_weights,
					XFLOAT op_significant_weight,    // TODO Put in const
					XFLOAT op_sum_weight,            // TODO Put in const
					unsigned long   coarse_trans,
					unsigned long   oversamples_trans,
					unsigned long   oversamples_orient,
					unsigned long   oversamples,
					bool  do_ignore_pdf_direction,
					XFLOAT *g_o_weights,
					XFLOAT *g_thr_wsum_prior_offsetx_class,
					XFLOAT *g_thr_wsum_prior_offsety_class,
					XFLOAT *g_thr_wsum_prior_offsetz_class,
					XFLOAT *g_thr_wsum_sigma2_offset,
					unsigned long *d_rot_idx,
					unsigned long *d_trans_idx,
					unsigned long *d_job_idx,
					unsigned long *d_job_num
				)
{
	// block id
	for (int bid=0; bid < grid_size; bid++) {

		XFLOAT s_o_weights[block_size];
		XFLOAT s_thr_wsum_sigma2_offset[block_size];;
		XFLOAT s_thr_wsum_prior_offsetx_class[block_size];
		XFLOAT s_thr_wsum_prior_offsety_class[block_size];
		XFLOAT s_thr_wsum_prior_offsetz_class[block_size];

		unsigned long pos = d_job_idx[bid];
		unsigned long job_size = d_job_num[bid];

		int pass_num = job_size/block_size + 1;

		for(int tid=0; tid<block_size; tid++) {
			s_o_weights[tid]                    	= (XFLOAT)0.0;
			s_thr_wsum_sigma2_offset[tid]       	= (XFLOAT)0.0;
			s_thr_wsum_prior_offsetx_class[tid] 	= (XFLOAT)0.0;
			s_thr_wsum_prior_offsety_class[tid] 	= (XFLOAT)0.0;
			if(DATA3D)
				s_thr_wsum_prior_offsety_class[tid] = (XFLOAT)0.0;
		}


		for (int pass = 0; pass < pass_num; pass++, pos+=block_size) // loop the available warps enough to complete all translations for this orientation
		{
			for(int tid=0; tid<block_size; tid++) {
				if ((pass*block_size+tid)<job_size) // if there is a translation that needs to be done still for this thread
				{
					// index of comparison
					long int iy = d_trans_idx[pos+tid];              // ...and its own trans...

					XFLOAT weight = g_i_weights[pos+tid];
					if( weight >= op_significant_weight ) //TODO Might be slow (divergent threads)
						weight /= op_sum_weight;
					else
						weight = (XFLOAT)0.0;

					s_o_weights[tid] += weight;
					s_thr_wsum_prior_offsetx_class[tid] += weight *          g_oo_otrans_x[iy];
					s_thr_wsum_prior_offsety_class[tid] += weight *          g_oo_otrans_y[iy];
					s_thr_wsum_sigma2_offset[tid]       += weight * g_myp_oo_otrans_x2y2z2[iy];
				}
			} 
		}

		for(int tid=1; tid<block_size; tid++)
		{
			s_o_weights[0]                    += s_o_weights[tid];
			s_thr_wsum_sigma2_offset[0]       += s_thr_wsum_sigma2_offset[tid];
			s_thr_wsum_prior_offsetx_class[0] += s_thr_wsum_prior_offsetx_class[tid];
			s_thr_wsum_prior_offsety_class[0] += s_thr_wsum_prior_offsety_class[tid];
			if(DATA3D)
				s_thr_wsum_prior_offsetz_class[0] += s_thr_wsum_prior_offsetz_class[tid];
		}
		g_o_weights[bid]			        = s_o_weights[0];
		g_thr_wsum_sigma2_offset[bid]       = s_thr_wsum_sigma2_offset[0];
		g_thr_wsum_prior_offsetx_class[bid] = s_thr_wsum_prior_offsetx_class[0];
		g_thr_wsum_prior_offsety_class[bid] = s_thr_wsum_prior_offsety_class[0];
		if(DATA3D)
			g_thr_wsum_prior_offsetz_class[bid] = s_thr_wsum_prior_offsetz_class[0];
	} // for bid
}

void exponentiate_weights_fine(
		XFLOAT *g_pdf_orientation,
		bool *g_pdf_orientation_zeros,
		XFLOAT *g_pdf_offset,
		bool *g_pdf_offset_zeros,
		XFLOAT *g_weights,
		XFLOAT min_diff2,
		unsigned long oversamples_orient,
		unsigned long oversamples_trans,
		unsigned long *d_rot_id,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		long int job_num);

void RNDnormalDitributionComplexWithPowerModulation2D(ACCCOMPLEX *Image, size_t xdim, XFLOAT *spectra);
void RNDnormalDitributionComplexWithPowerModulation3D(ACCCOMPLEX *Image, size_t xdim, size_t ydim, XFLOAT *spectra);

void softMaskBackgroundValue(	int      block_dim,
                                int      block_size,
                                XFLOAT  *vol,
								long int vol_size,
								long int xdim,
								long int ydim,
								long int zdim,
								long int xinit,
								long int yinit,
								long int zinit,
								XFLOAT   radius,
								XFLOAT   radius_p,
								XFLOAT   cosine_width,
								XFLOAT  *g_sum,
								XFLOAT  *g_sum_bg);

void cosineFilter(	int      block_dim,
                    int      block_size,
					XFLOAT *vol,
					long int vol_size,
					long int xdim,
					long int ydim,
					long int zdim,
					long int xinit,
					long int yinit,
					long int zinit,
					bool do_Mnoise,
					XFLOAT *noise,
					XFLOAT radius,
					XFLOAT radius_p,
					XFLOAT cosine_width,
					XFLOAT sum_bg_total);
								
//----------------------------------------------------------------------------

template <typename T>
void sycl_gpu_translate2D(T *g_image_in,
					T		*g_image_out,
					size_t   image_size,
					int      xdim,
					int      ydim, //not used
					int      dx,
					int      dy,
					virtualSYCL	*devAcc)
{
	assert(image_size <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		cgh.parallel_for(sycl::range<1>(image_size), [=](sycl::id<1> pixel)
		{
			int x = pixel % xdim;
			int y = (pixel-x) / xdim;

			int xp = x + dx;
			int yp = y + dy;

			if( yp>=0 && xp>=0 && yp<ydim && xp<xdim)
			{
				int new_pixel = yp*xdim + xp;
				if(new_pixel>=0 && new_pixel<image_size) // if displacement is negative, new_pixel could be less than 0
					g_image_out[new_pixel] = g_image_in[pixel];
			}

		});
	}).wait_and_throw();
}

template <typename T>
void sycl_gpu_translate3D(T *g_image_in,
					T		*g_image_out,
					size_t   image_size,
					int      xdim,
					int      ydim,
					int      zdim, //not used
					int      dx,
					int      dy,
					int      dz,
					virtualSYCL	*devAcc)
{
	assert(image_size <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		cgh.parallel_for(sycl::range<1>(image_size), [=](sycl::id<1> voxel)
		{
			int xydim = xdim*ydim;

			int z =  voxel / xydim;
			int zp = z + dz;

			int xy = voxel % xydim;
			int y =  xy / xdim;
			int yp = y + dy;

			int x =  xy % xdim;
			int xp = x + dx;

			if( zp>=0 && yp>=0 && xp>=0 && zp<zdim && yp<ydim && xp<xdim)
			{
				int new_voxel = zp*xydim +  yp*xdim + xp;
				if(new_voxel>=0 && new_voxel<image_size) // if displacement is negative, new_pixel could be less than 0
					g_image_out[new_voxel] = g_image_in[voxel];
			}
		});
	}).wait_and_throw();
}

//----------------------------------------------------------------------------
/*
 * Multiplies scalar array A by a scalar S
 *
 *  OUT[i] = A[i]*S
 */
template <typename T>
void sycl_gpu_kernel_multi( T   *A,
            T   *OUT,
            T    S,
            size_t   image_size,
			virtualSYCL *devAcc)
{
	assert(image_size <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		cgh.parallel_for(sycl::range<1>(image_size), [=](sycl::id<1> idx)
		{
			OUT[idx] = A[idx]*S;
		});
	}).wait_and_throw();
}
/*
 * In place multiplies scalar array A by a scalar S
 *
 *  A[i] = A[i]*S
 */
template <typename T>
void sycl_gpu_kernel_multi( T   *A,
			T    S,
			size_t     image_size,
			virtualSYCL *devAcc)
{
	assert(image_size <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		cgh.parallel_for(sycl::range<1>(image_size), [=](sycl::id<1> idx)
		{
			A[idx] *= S;
		});
	}).wait_and_throw();
}
/*
 * Multiplies scalar array A by scalar array B and a scalar S, pixel-by-pixel
 *
 *  OUT[i] = A[i]*B[i]*S
 */
template <typename T>
void sycl_gpu_kernel_multi( T *A,
			T *B,
			T *OUT,
			T  S,
			size_t   image_size,
			virtualSYCL *devAcc)
{
	assert(image_size <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		cgh.parallel_for(sycl::range<1>(image_size), [=](sycl::id<1> idx)
		{
			OUT[idx] = A[idx]*B[idx]*S;
		});
	}).wait_and_throw();
}

/*
 * In place add scalar S to scalar array A
 *
 *  A[i] = A[i]*S
 */
template <typename T>
void sycl_gpu_kernel_add(
    T *A,
    T S,
    size_t size,
	virtualSYCL *devAcc
)
{
	assert(size <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);

	dev->syclSubmit([&](sycl::handler &cgh)
	{
		cgh.parallel_for(sycl::range<1>(size), [=](sycl::id<1> idx)
		{
			A[idx] += S;
		});
	}).wait_and_throw();
}

template<bool do_highpass>
inline void kernel_frequencyPass( int grid_size, int block_size,
					ACCCOMPLEX *A,
					long int     ori_size,
					size_t       Xdim,
					size_t       Ydim,
					size_t       Zdim,
					XFLOAT       edge_low,
					XFLOAT       edge_width,
					XFLOAT       edge_high,
					XFLOAT       angpix,
					size_t       image_size)
{
	// TODO - why not a single loop over image_size pixels?
	for(int blk=0; blk<grid_size; blk++) {
		for(int tid=0; tid<block_size; tid++) {
			size_t texel = (size_t)tid + (size_t)blk*(size_t)block_size;

			int z = texel / (Xdim*Ydim);
			int xy = (texel - z*Xdim*Ydim);
			int y = xy / Xdim;

			int xp = xy - y*Xdim;

			int zp = ( z<Xdim ? z : z-Zdim );
			int yp = ( y<Xdim ? y : y-Ydim );

			int r2 = xp*xp + yp*yp + zp*zp;

			RFLOAT res;
			if(texel<image_size)
			{
				res = sqrt((RFLOAT)r2)/(RFLOAT)ori_size;

				if(do_highpass) //highpass
				{
					if (res < edge_low) //highpass => lows are dead
					{
						A[texel].x = 0.;
						A[texel].y = 0.;
					}
					else if (res < edge_high) //highpass => medium lows are almost dead
					{
						XFLOAT mul = 0.5 - 0.5 * std::cos( PI * (res-edge_low)/edge_width);
						A[texel].x *= mul;
						A[texel].y *= mul;
					}
				}
				else //lowpass
				{
					if (res > edge_high) //lowpass => highs are dead
					{
						A[texel].x = 0.;
						A[texel].y = 0.;
					}
					else if (res > edge_low) //lowpass => medium highs are almost dead
					{
						XFLOAT mul = 0.5 + 0.5 * std::cos( PI * (res-edge_low)/edge_width);
						A[texel].x *= mul;
						A[texel].y *= mul;
					}
				}
			}
		} // tid
	} // blk
}

template<bool DATA3D>
inline void powerClass(int          gridSize,
				ACCCOMPLEX   *g_image,
				XFLOAT       *g_spectrum,
				size_t       image_size,
				size_t       spectrum_size,
				int          xdim,
				int          ydim,
				int          zdim,
				int          res_limit,
				XFLOAT      *g_highres_Xi2)
{
	for(int bid=0; bid<gridSize; bid++)
	{
		XFLOAT normFaux;

		int x,y,xy,d;
		int xydim = xdim*ydim;
		bool coords_in_range(true);

		XFLOAT s_highres_Xi2[POWERCLASS_BLOCK_SIZE];
		for(int tid=0; tid<POWERCLASS_BLOCK_SIZE; tid++)
			s_highres_Xi2[tid] = (XFLOAT)0.;

		for(int tid=0; tid<POWERCLASS_BLOCK_SIZE; tid++){
			size_t voxel=(size_t)tid + (size_t)bid*(size_t)POWERCLASS_BLOCK_SIZE;
			if(voxel<image_size)
			{
				if(DATA3D)
				{
					int z =  voxel / xydim;
					xy = voxel % xydim;
					y =  xy / xdim;
					x =  xy % xdim;

					y = ((y<xdim) ? y : y-ydim);
					z = ((z<xdim) ? z : z-zdim);

					d  = (x*x + y*y + z*z);
					coords_in_range = !(x==0 && y<0.f && z<0.f);
				}
				else
				{
					x = voxel % xdim;
					y = (voxel-x) / (xdim);

					y = ((y<xdim) ? y : y-ydim);
					d  = (x*x + y*y);
					coords_in_range = !(x==0 && y<0.f);
				}

				size_t ires = static_cast<size_t>(std::sqrt(static_cast<XFLOAT>(d)) + 0.5);
				if((ires<spectrum_size) && coords_in_range)
				{
					normFaux = g_image[voxel].x*g_image[voxel].x + g_image[voxel].y*g_image[voxel].y;
					g_spectrum[ires] += normFaux;
					if(ires>=res_limit)
						s_highres_Xi2[tid] += normFaux;
				}
			}
		}

		for(int tid=1; tid<POWERCLASS_BLOCK_SIZE; tid++)
			s_highres_Xi2[0] += s_highres_Xi2[tid];

		g_highres_Xi2[0] += s_highres_Xi2[0];
	}
}

inline std::pair<XFLOAT, XFLOAT> sycl_sincos(XFLOAT val)
{
#ifdef ACC_DOUBLE_PRECISION
	return std::make_pair(sycl::sin(val), sycl::cos(val));
#else
	return std::make_pair(sycl::native::sin(val), sycl::native::cos(val));
#endif
}

inline void translatePixel(
		int x,
		int y,
		XFLOAT tx,
		XFLOAT ty,
		XFLOAT real,
		XFLOAT imag,
		XFLOAT &tReal,
		XFLOAT &tImag)
{
	XFLOAT v = x * tx + y * ty;
#ifdef ACC_DOUBLE_PRECISION
	XFLOAT s = sycl::sin(v);
	XFLOAT c = sycl::cos(v);
#else
	XFLOAT s = sycl::native::sin(v);
	XFLOAT c = sycl::native::cos(v);
#endif

	tReal = c * real - s * imag;
	tImag = c * imag + s * real;
}

inline void translatePixel(
		int x,
		int y,
		int z,
		XFLOAT tx,
		XFLOAT ty,
		XFLOAT tz,
		XFLOAT real,
		XFLOAT imag,
		XFLOAT &tReal,
		XFLOAT &tImag)
{
	XFLOAT v = x * tx + y * ty + z * tz;
#ifdef ACC_DOUBLE_PRECISION
	XFLOAT s = sycl::sin(v);
	XFLOAT c = sycl::cos(v);
#else
	XFLOAT s = sycl::native::sin(v);
	XFLOAT c = sycl::native::cos(v);
#endif

	tReal = c * real - s * imag;
	tImag = c * imag + s * real;
}

// sincos lookup table optimization. Function translatePixel calls 
// sincos(x*tx + y*ty). We precompute 2D lookup tables for x and y directions. 
// The first dimension is x or y pixel index, and the second dimension is x or y
// translation index. Since sin(a+B) = sin(A) * cos(B) + cos(A) * sin(B), and 
// cos(A+B) = cos(A) * cos(B) - sin(A) * sin(B), we can use lookup table to 
// compute sin(x*tx + y*ty) and cos(x*tx + y*ty). 
inline void  computeSincosLookupTable2D(unsigned long  trans_num,
                                        XFLOAT  *trans_x,
										XFLOAT  *trans_y,										
										int      xSize,
										int      ySize,
										XFLOAT  *sin_x,
										XFLOAT  *cos_x,
	                                    XFLOAT  *sin_y,
	                                    XFLOAT  *cos_y)
{
	for(unsigned long i=0; i<trans_num; i++) {
		XFLOAT tx = trans_x[i];
		XFLOAT ty = trans_y[i];

		for(int x=0; x<xSize; x++) {
			unsigned long index = i * xSize + x;
			XFLOAT v = x * tx;
#ifdef ACC_DOUBLE_PRECISION
			sin_x[index] = sycl::sin(v);
			cos_x[index] = sycl::cos(v);
#else
			sin_x[index] = sycl::native::sin(v);
			cos_x[index] = sycl::native::cos(v);
#endif
		}
		
		for(int y=0; y<ySize; y++) {
			unsigned long index = i * ySize + y;
			XFLOAT v = y * ty;
#ifdef ACC_DOUBLE_PRECISION
			sin_y[index] = sycl::sin(v);
			cos_y[index] = sycl::cos(v);
#else
			sin_y[index] = sycl::native::sin(v);
			cos_y[index] = sycl::native::cos(v);
#endif
		}
	}
}
				
inline void  computeSincosLookupTable3D(unsigned long  trans_num,
                                        XFLOAT  *trans_x,
										XFLOAT  *trans_y,
										XFLOAT  *trans_z,										
										int      xSize,
										int      ySize,
										int      zSize,
										XFLOAT  *sin_x,
										XFLOAT  *cos_x,
	                                    XFLOAT  *sin_y,
	                                    XFLOAT  *cos_y,
	                                    XFLOAT  *sin_z,
	                                    XFLOAT  *cos_z)
{	                                    
	for(unsigned long i=0; i<trans_num; i++) {
		XFLOAT tx = trans_x[i];
		XFLOAT ty = trans_y[i];
		XFLOAT tz = trans_z[i];

		for(int x=0; x<xSize; x++) {
			unsigned long index = i * xSize + x;
			XFLOAT v = x * tx;
#ifdef ACC_DOUBLE_PRECISION
			sin_x[index] = sycl::sin(v);
			cos_x[index] = sycl::cos(v);
#else
			sin_x[index] = sycl::native::sin(v);
			cos_x[index] = sycl::native::cos(v);
#endif
		}
		
		for(int y=0; y<ySize; y++) {
			unsigned long index = i * ySize + y;
			XFLOAT v = y * ty;
#ifdef ACC_DOUBLE_PRECISION
			sin_y[index] = sycl::sin(v);
			cos_y[index] = sycl::cos(v);
#else
			sin_y[index] = sycl::native::sin(v);
			cos_y[index] = sycl::native::cos(v);
#endif
		}
		
		for(int z=0; z<zSize; z++) {
			unsigned long index = i * zSize + z;
			XFLOAT v = z * tz;
#ifdef ACC_DOUBLE_PRECISION
			sin_z[index] = sycl::sin(v);
			cos_z[index] = sycl::cos(v);
#else
			sin_z[index] = sycl::native::sin(v);
			cos_z[index] = sycl::native::cos(v);
#endif
		}
	}
}

template<bool invert>
void sycl_kernel_make_eulers_2D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *eulers,
		unsigned long orientation_num);

template<bool invert,bool doL, bool doR>
void sycl_kernel_make_eulers_3D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *betas,
		XFLOAT *gammas,
		XFLOAT *eulers,
		unsigned long orientation_num,
		XFLOAT *L,
		XFLOAT *R);

template< typename T>
size_t findThresholdIdxInCumulativeSum(T *data, const size_t count, T threshold, virtualSYCL *devAcc)
{
	assert(count <= std::numeric_limits<int>::max());
	sycl::queue *Q = dynamic_cast<devSYCL*>(devAcc)->getQueue();
	size_t dist;
 #ifdef USE_MORE_ONEDPL
	auto iter = dpl::find_if(dynamic_cast<devSYCL*>(devAcc)->getDevicePolicy(), data, data+count, [=](T &v){return (v>threshold);});
	dist = dpl::distance(data, iter);
 #else
	dist = 0;
	{
		using namespace sycl;
		auto buf = buffer{&dist, range{1}};

		Q->submit([&](handler &cgh)
		{
			auto loc = accessor{buf, cgh, write_only};
			cgh.parallel_for(range<1>(count-1), [=](id<1> idx)
			{
				if (data[idx] <= threshold && data[idx+1] > threshold)
					loc[0] = idx+1;
			});
		}).wait_and_throw();
	}
 #endif
	return dist;
}

static void sycl_kernel_allweights_to_mweights(unsigned long *d_iorient, XFLOAT *d_allweights, XFLOAT *d_mweights, unsigned long orientation_num, unsigned long translation_num, int block_size, virtualSYCL *devAcc)
{
	assert(orientation_num*translation_num <= std::numeric_limits<int>::max());
	auto dev = dynamic_cast<devSYCL*>(devAcc);
	sycl::event prev;
	bool isAsync = dev->isAsyncQueue();
	if (isAsync)
		prev = dev->getLastEvent();

	auto ev = dev->syclSubmit([&](sycl::handler &cgh)
	{
		using namespace sycl;
		if (isAsync)
			cgh.depends_on(prev);

		cgh.parallel_for(range<1>(orientation_num*translation_num), [=](id<1> idx)
		{
			d_mweights[d_iorient[idx/translation_num] * translation_num + idx%translation_num]
				=
					d_allweights[idx/translation_num * translation_num + idx%translation_num];
		});
	});
	dev->pushEvent(ev);
}

} // end of namespace syclGpuKernels

#endif /* HELPER_KERNELS_H_ */
