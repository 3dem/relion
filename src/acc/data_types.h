#ifndef ACC_DATA_TYPES_H_
#define ACC_DATA_TYPES_H_

#include "src/acc/acc_ptr.h"
#include "src/multidim_array.h"

namespace AccDataTypes
{

template <typename T, int AccT = ACC_CPU>
class Image : public AccPtr<T,AccT>
{
private:
	int x,y,z;
	bool fourier; //Is this a Fourier space data array

public:

	/*======================================================
				         CONSTRUCTORS
	======================================================*/

	Image(CudaCustomAllocator *allocator):
		AccPtr<T,AccT>(allocator),
		x(0), y(0), z(0), fourier(false)
	{}

	Image(cudaStream_t stream, CudaCustomAllocator *allocator):
		AccPtr<T,AccT>(stream, allocator),
		x(0), y(0), z(0), fourier(false)
	{}

	Image(int xdim, int ydim, CudaCustomAllocator *allocator):
		AccPtr<T,AccT>(xdim*ydim, allocator),
		x(xdim), y(ydim), z(1), fourier(false)
	{}

	Image(int xdim, int ydim, int zdim, CudaCustomAllocator *allocator):
		AccPtr<T,AccT>(xdim*ydim*zdim, allocator),
		x(xdim), y(ydim), z(zdim), fourier(false)
	{}

	Image(int xdim, int ydim, cudaStream_t stream, CudaCustomAllocator *allocator):
		 AccPtr<T,AccT>(xdim*ydim, stream, allocator),
		 x(xdim), y(ydim), z(1), fourier(false)
	{}

	Image(int xdim, int ydim, int zdim, cudaStream_t stream, CudaCustomAllocator *allocator):
		 AccPtr<T,AccT>(xdim*ydim*zdim, stream, allocator),
		 x(xdim), y(ydim), z(zdim), fourier(false)
	{}

	template <typename T1>
	Image(MultidimArray<T1> img, CudaCustomAllocator *allocator):
		AccPtr<T,AccT>(img.nzyxdim, allocator),
		x(img.xdim), y(img.ydim), z(img.zdim), fourier(false)
	{}

	template <typename T1>
	Image(MultidimArray<T1> img, cudaStream_t stream, CudaCustomAllocator *allocator):
		 AccPtr<T,AccT>(img.nzyxdim, stream, allocator),
		x(img.xdim), y(img.ydim), z(img.zdim), fourier(false)
	{}

	Image(int box_dim, bool is_fourier, bool is3D, CudaCustomAllocator *allocator)
	{
		setSize(box_dim, is_fourier, is3D);
		AccPtr<T,AccT>(x*y*z, allocator);
	}

	Image(int box_dim, bool is_fourier, bool is3D, cudaStream_t stream, CudaCustomAllocator *allocator)
	{
		setSize(box_dim, is_fourier, is3D);
		AccPtr<T,AccT>(x*y*z, stream, allocator);
	}

	/*======================================================
				           METHODS
	======================================================*/

	int getx() { return x; }
	int gety() { return y; }
	int getz() { return z; }

	int getxy()
	{
		return x*y;
	}

	int getxyz()
	{
		return AccPtr<T,AccT>::getSize();
	}

	bool is3D() { return z > 1; }

	void setSize(int box_dim, bool is_fourier, bool is3D)
	{
		fourier = is_fourier;
		if (is_fourier)
		{
			x = box_dim/2+1;
			y = box_dim;
			if (is3D)
				z = box_dim;
			else
				z = 1;
		}
		AccPtr<T,AccT>::setSize(x*y*z);
	}

	void setSize(int xdim, int ydim)
	{
		x = xdim;
		y = ydim;

		AccPtr<T,AccT>::setSize(x*y);
	}

	void setSize(int xdim, int ydim, int zdim)
	{
		x = xdim;
		y = ydim;
		z = zdim;

		AccPtr<T,AccT>::setSize(x*y*z);
	}

	template <typename T1>
	void setSize(MultidimArray<T1> img)
	{
		x = img.xdim;
		y = img.xdim;
		z = img.xdim;

		AccPtr<T,AccT>::setSize(x*y*z);
	}

	template <typename T1>
	void setHost(MultidimArray<T1> img)
	{
		if (img.xdim != x || img.ydim != y || img.zdim != z)
		{
			if (img.nzyxdim > AccPtr<T,AccT>::getSize())
			{
				AccPtr<T,AccT>::freeIfSet();
				setSize(img);
				AccPtr<T,AccT>::hostAlloc();
			}
			else
				setSize(img);
		}

		T *ptr = AccPtr<T,AccT>::getHostPtr();

		if (sizeof(T) == sizeof(T1))
			memcpy(ptr, img.data, sizeof(T)*img.nzyxdim);
		else
			for (int i = 0; i < img.nzyxdim; i++)
				ptr[i] = (T) img.data[i];
	}
};

}

#endif
