#ifndef ACC_BACKPROJECTOR_H_
#define ACC_BACKPROJECTOR_H_

#ifdef _CUDA_ENABLED
#  include <cuda_runtime.h>
using deviceStream_t = cudaStream_t;
#elif _HIP_ENABLED
#  include <hip/hip_runtime.h>
using deviceStream_t = hipStream_t;
#elif _SYCL_ENABLED
#  include "src/acc/sycl/sycl_virtual_dev.h"
using deviceStream_t = virtualSYCL*;
#else
using deviceStream_t = float;
#endif
#include "src/complex.h"
#include "src/acc/settings.h"
#include "src/acc/acc_ptr.h"

#ifdef ALTCPU
#  include <tbb/spin_mutex.h>
#endif

class AccBackprojector
{

public:
	int mdlX, mdlY, mdlZ,
	    mdlInitY, mdlInitZ,
	    maxR, maxR2;
	XFLOAT padding_factor;
	size_t mdlXYZ;

#ifdef ALTCPU
	tbb::spin_mutex *mutexes;
#endif

	size_t allocaton_size;
	size_t voxelCount;

	XFLOAT *d_mdlReal, *d_mdlImag, *d_mdlWeight;

	deviceStream_t stream;

public:

	AccBackprojector():
				mdlX(0), mdlY(0), mdlZ(0), mdlXYZ(0),
				mdlInitY(0), mdlInitZ(0),
				maxR(0), maxR2(0),
				padding_factor(0),
				allocaton_size(0), voxelCount(0),
				d_mdlReal(NULL), d_mdlImag(NULL), d_mdlWeight(NULL),
				stream(0)
#ifdef ALTCPU
				, mutexes(0)
#endif
	{}

	size_t setMdlDim(
#ifdef _SYCL_ENABLED
			deviceStream_t dev,
#endif
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, XFLOAT paddingFactor);

	void initMdl();

	void backproject(
			XFLOAT *d_imgs_nomask_real,
			XFLOAT *d_imgs_nomask_imag,
			XFLOAT *trans_x,
			XFLOAT *trans_y,
			XFLOAT *trans_z,
			XFLOAT* d_weights,
			XFLOAT* d_Minvsigma2s,
			XFLOAT* d_ctfs,
			unsigned long translation_num,
			XFLOAT significant_weight,
			XFLOAT weight_norm,
			XFLOAT *d_eulers,
			int imgX,
			int imgY,
			int imgZ,
			unsigned long imageCount,
			bool data_is_3D,
			deviceStream_t optStream);

	void getMdlData(XFLOAT *real, XFLOAT *imag, XFLOAT * weights);
	void getMdlDataPtrs(XFLOAT *& real, XFLOAT *& imag, XFLOAT *& weights);

	void setStream(deviceStream_t s) { stream = s; }
	deviceStream_t getStream() { return stream; }

	void clear();

	~AccBackprojector();
};

#endif
