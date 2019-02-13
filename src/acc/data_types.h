#ifndef ACC_DATA_TYPES_H_
#define ACC_DATA_TYPES_H_

#include "src/acc/acc_ptr.h"
#include "src/multidim_array.h"

namespace AccDataTypes
{

template <typename T>
class Image : public AccPtr<T>
{
private:
	int x,y,z;
	bool fourier; //Is this a Fourier space data array

public:

	/*======================================================
				         CONSTRUCTORS
	======================================================*/

	Image(AccPtrFactory &f):
		AccPtr<T>(f.make<T>()),
		x(0), y(0), z(0), fourier(false)
	{}

	Image(int xdim, AccPtrFactory &f):
		AccPtr<T>(f.make<T>(xdim)),
		x(xdim), y(1), z(1), fourier(false)
	{}

	Image(int xdim, int ydim, AccPtrFactory &f):
		AccPtr<T>(f.make<T>(xdim*ydim)),
		x(xdim), y(ydim), z(1), fourier(false)
	{}

	Image(int xdim, int ydim, int zdim, AccPtrFactory &f):
		AccPtr<T>(f.make<T>(xdim*ydim*zdim)),
		x(xdim), y(ydim), z(zdim), fourier(false)
	{}

	template <typename T1>
	Image(MultidimArray<T1> img, AccPtrFactory &f):
		AccPtr<T>(f.make<T>(img.nzyxdim)),
		x(img.xdim), y(img.ydim), z(img.zdim), fourier(false)
	{}

	Image(int box_dim, bool is_fourier, bool is3D, AccPtrFactory &f)
	{
		setSize(box_dim, is_fourier, is3D);
		AccPtr<T>(f.make<T>(x*y*z));
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
		return AccPtr<T>::getSize();
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
		AccPtr<T>::setSize(x*y*z);
	}

	void setSize(int xdim)
	{
		x = xdim;
		y = 1;
		z = 1;

		AccPtr<T>::setSize(x);
	}

	void setSize(int xdim, int ydim)
	{
		x = xdim;
		y = ydim;
		z = 1;

		AccPtr<T>::setSize(x*y);
	}

	void setSize(int xdim, int ydim, int zdim)
	{
		x = xdim;
		y = ydim;
		z = zdim;

		AccPtr<T>::setSize(x*y*z);
	}

	template <typename T1>
	void setSize(MultidimArray<T1> img)
	{
		x = img.xdim;
		y = img.xdim;
		z = img.xdim;

		AccPtr<T>::setSize(x*y*z);
	}

	template <typename T1>
	void setHost(MultidimArray<T1> &img)
	{
		if (img.xdim != x || img.ydim != y || img.zdim != z)
		{
			if (img.nzyxdim > AccPtr<T>::getSize())
			{
				AccPtr<T>::freeIfSet();
				setSize(img);
				AccPtr<T>::hostAlloc();
			}
			else
				setSize(img);
		}

		if (AccPtr<T>::getHostPtr() == NULL)
			AccPtr<T>::hostAlloc();

		T *ptr = AccPtr<T>::getHostPtr();

		if (sizeof(T) == sizeof(T1))
			memcpy(ptr, img.data, sizeof(T)*img.nzyxdim);
		else
			for (unsigned long i = 0; i < img.nzyxdim; i++)
				ptr[i] = (T) img.data[i];
	}

	template <typename T1>
	void getHost(MultidimArray<T1> &img)
	{

		if(img.nzyxdim!=AccPtr<T>::getSize())
		{
			if(img.nzyxdim==0)
				img.resize(z,y,x);
			else
				CRITICAL("Trying to fill host-array with data from an array with different size!")
		}
		T *ptr = AccPtr<T>::getHostPtr();

		if (sizeof(T) == sizeof(T1))
			memcpy(img.data, ptr, sizeof(T)*img.nzyxdim);
		else
			for (unsigned long i = 0; i < img.nzyxdim; i++)
				img.data[i] = (T1) ptr[i];
	}
};

}

#endif
