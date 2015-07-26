#ifndef CUDA_BACKPROJECTOR_H_
#define CUDA_BACKPROJECTOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"
#include <cuda_runtime.h>

#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_FETCH_COUNT 4

class CudaBackprojector
{
	int mdlX, mdlY, mdlZ, mdlXYZ,
	    mdlInitY, mdlInitZ,
	    maxR2,
	    padding_factor;

	int *h_voxelX, *h_voxelY, *h_voxelZ;

	int *d_voxelX, *d_voxelY, *d_voxelZ;
	FLOAT *d_mdlReal, *d_mdlImag, *d_mdlWeight;

	unsigned long voxelCount;

	cudaStream_t stream;

public:

	CudaBackprojector():
			mdlX(0), mdlY(0), mdlZ(0), mdlXYZ(0),
			mdlInitY(0), mdlInitZ(0),
			maxR2(0),
			padding_factor(0),
			h_voxelX(0), h_voxelY(0), h_voxelZ(0),
			d_voxelX(0), d_voxelY(0), d_voxelZ(0),
			voxelCount(0),
			d_mdlReal(0), d_mdlImag(0), d_mdlWeight(0),
			stream(0)
	{}

	CudaBackprojector(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int padding_factor):
			mdlX(xdim), mdlY(ydim), mdlZ(zdim),
			mdlXYZ(xdim*ydim*zdim),
			mdlInitY(inity), mdlInitZ(initz),
			maxR2(max_r*max_r),
			padding_factor(padding_factor),
			h_voxelX(0), h_voxelY(0), h_voxelZ(0),
			d_voxelX(0), d_voxelY(0), d_voxelZ(0),
			voxelCount(0),
			d_mdlReal(0), d_mdlImag(0), d_mdlWeight(0),
			stream(0)
	{}

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

	cudaStream_t getStream() { return stream; }

	~CudaBackprojector();

};

#endif
