#ifndef CUDA_BACKPROJECTOR_H_
#define CUDA_BACKPROJECTOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"

#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_FETCH_COUNT 4

class Cuda3DBackprojector
{
public:
	int mdlX, mdlY, mdlZ, mdlXYZ,
	    mdlInitY, mdlInitZ,
	    maxR2,
	    padding_factor;

	void *voxelX, *voxelY, *voxelZ; //CudaGlobalPtr<int>
	unsigned long voxelCount;

	void *mdlReal, *mdlImag, *mdlWeight; //CudaGlobalPtr<FLOAT>

	Cuda3DBackprojector():
			mdlX(0), mdlY(0), mdlZ(0), mdlXYZ(0),
			mdlInitY(0), mdlInitZ(0),
			maxR2(0),
			padding_factor(0),
			voxelX(0), voxelY(0), voxelZ(0),
			voxelCount(0),
			mdlReal(0), mdlImag(0), mdlWeight(0)
	{};

	Cuda3DBackprojector(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int padding_factor):
			mdlX(xdim), mdlY(ydim), mdlZ(zdim),
			mdlXYZ(xdim*ydim*zdim),
			mdlInitY(inity), mdlInitZ(initz),
			maxR2(max_r*max_r),
			padding_factor(padding_factor),
			voxelX(0), voxelY(0), voxelZ(0),
			voxelCount(0),
			mdlReal(0), mdlImag(0), mdlWeight(0)
	{};

	void setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int paddingFactor)
	{
		mdlX = xdim;
		mdlY = ydim;
		mdlZ = zdim;
		mdlXYZ = xdim*ydim*zdim;
		mdlInitY = inity;
		mdlInitZ = initz;
		maxR2 = max_r*max_r;
		padding_factor = paddingFactor;
	}

	void initMdl();

	void backproject(
			FLOAT *d_real,
			FLOAT *d_imag,
			FLOAT *d_weight,
			FLOAT *d_eulers,
			int imgX,
			int imgY,
			unsigned long imageCount);

	void getMdlData(FLOAT *real, FLOAT *imag, FLOAT * weights);
	void getMdlData(Complex *data, double * weights);

	~Cuda3DBackprojector();

};

#endif
