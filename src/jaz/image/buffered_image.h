#ifndef JAZ_VECTOR_IMAGE_H
#define JAZ_VECTOR_IMAGE_H

#include "raw_image.h"


template <typename T>
class BufferedImage : public RawImage<T>
{
	public:
		
		template <typename T2>
		static BufferedImage<T> convert(const BufferedImage<T2>& img);
		
		BufferedImage();
		BufferedImage(size_t xdim, size_t ydim = 1, size_t zdim = 1);
		BufferedImage(const BufferedImage& vi);
		BufferedImage(const RawImage<T>& vi);
		BufferedImage(const Image<T>& vi);
		BufferedImage(std::string filename);
		
		
			std::vector<T> dataVec;
		
			
		void resize(size_t xdim, size_t ydim, size_t zdim = 1);
		void resize(const RawImage<T>& example);
		
		void read(std::string fn, int slice = -1);
		
		template <typename T2>
		void copyDataAndSizeFrom(const gravis::tImage<T2>& img);
		
		template <typename T2>
		void copyDataAndSizeFrom(const Image<T2>& img);
		
		RawImage<T> getRef();
		
		BufferedImage& operator = (const BufferedImage& other);
};


template <typename T>
template <typename T2>
BufferedImage<T> BufferedImage<T>::convert(const BufferedImage<T2>& img)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	const size_t s = w * h * d;
	
	BufferedImage<T> out(w,h,d);
	
	for (int i = 0; i < s; i++)
	{
		out[i] = (T) img[i];
	}
	
	return out;
}


template <class T>
BufferedImage<T>& BufferedImage<T>::operator = (const BufferedImage<T>& other)
{
	RawImage<T>::operator=(other);
	
	dataVec = other.dataVec;
	RawImage<T>::data = &(dataVec[0]);

	return *this;
}

template <class T>
BufferedImage<T>::BufferedImage()
	: RawImage<T>()
{}

template <class T>
BufferedImage<T>::BufferedImage(size_t xdim, size_t ydim, size_t zdim)
	:   RawImage<T>(xdim, ydim, zdim, 0)
{
	dataVec.resize(xdim*ydim*zdim);
	RawImage<T>::data = &(dataVec[0]);
}

template <class T>
BufferedImage<T>::BufferedImage(const BufferedImage<T>& vi)
	:   RawImage<T>(vi),
	  dataVec(vi.dataVec)
{
	RawImage<T>::data = &(dataVec[0]);
}

template <class T>
BufferedImage<T>::BufferedImage(const RawImage<T>& vi)
	:   RawImage<T>(vi)
{
	const long int s = this->xdim * this->ydim * this->zdim;

	dataVec.resize(s);

	for (int i = 0; i < s; i++)
	{
		dataVec[i] = vi[i];
	}

	RawImage<T>::data = &(dataVec[0]);
}

template <class T>
BufferedImage<T>::BufferedImage(const Image<T>& vi)
	:   RawImage<T>(
			vi.data.xdim,
			vi.data.ydim,
			vi.data.zdim * vi.data.ndim,
			 0)
{
	const long int s = this->xdim * this->ydim * this->zdim;

	dataVec.resize(s);

	for (int i = 0; i < s; i++)
	{
		dataVec[i] = vi.data.data[i];
	}

	RawImage<T>::data = &(dataVec[0]);
}

template<typename T>
BufferedImage<T>::BufferedImage(std::string filename)
{
	read(filename);
}

template <class T>
void BufferedImage<T>::resize(size_t xdim, size_t ydim, size_t zdim)
{
	this->xdim = xdim;
	this->ydim = ydim;
	this->zdim = zdim;
	
	dataVec.resize(xdim*ydim*zdim);
	RawImage<T>::data = &(dataVec[0]);
}

template <class T>
void BufferedImage<T>::resize(const RawImage<T>& example)
{
	this->xdim = example.xdim;
	this->ydim = example.ydim;
	this->zdim = example.zdim;
	
	dataVec.resize(this->xdim * this->ydim * this->zdim);
	this->data = &(dataVec[0]);
}

template <class T>
void BufferedImage<T>::read(std::string fn, int slice)
{
	std::string::size_type dot = fn.find_last_of('.');
	
	if (dot == std::string::npos)
	{
		REPORT_ERROR_STR("Image<T>::read: filename has no ending (" << fn << ").");
	}
	
	std::string end = fn.substr(dot+1);
	
	if (end == "png")
	{
		gravis::tImage<T> img;
		img.read(fn);
		copyDataAndSizeFrom(img);
	}
	else
	{
		Image<T> img;
		img.read(fn, true, slice);
		copyDataAndSizeFrom(img);
	}	
}


template <class T>
template <typename T2>
inline void BufferedImage<T>::copyDataAndSizeFrom(const gravis::tImage<T2>& img)
{
	this->xdim = img.cols();
	this->ydim = img.rows();
	this->zdim = 1;
	
	dataVec.resize(this->xdim * this->ydim * this->zdim);
	RawImage<T>::data = &(dataVec[0]);
	
	for (size_t y = 0; y < this->ydim; y++)
	for (size_t x = 0; x < this->xdim; x++)
	{
		RawImage<T>::data[y * this->xdim + x] = (T)img(x,y);
	}
}

template <class T>
template <typename T2>
inline void BufferedImage<T>::copyDataAndSizeFrom(const Image<T2>& img)
{
	this->xdim = img.data.xdim;
	this->ydim = img.data.ydim;
	
	if (img.data.ndim == 1)
	{
		this->zdim = img.data.zdim;
		
		dataVec.resize(this->xdim * this->ydim * this->zdim);
		RawImage<T>::data = &(dataVec[0]);
		
		for (size_t z = 0; z < this->zdim; z++)
		for (size_t y = 0; y < this->ydim; y++)
		for (size_t x = 0; x < this->xdim; x++)
		{
			T value = (T)(DIRECT_A3D_ELEM(img(), z, y, x));
			RawImage<T>::data[(z * this->ydim + y) * this->xdim + x] = value;
		}
	}
	else
	{
		this->zdim = img.data.ndim;
		
		dataVec.resize(this->xdim * this->ydim * this->zdim);
		RawImage<T>::data = &(dataVec[0]);
		
		for (size_t z = 0; z < this->zdim; z++)
		for (size_t y = 0; y < this->ydim; y++)
		for (size_t x = 0; x < this->xdim; x++)
		{
			T value = (T)(DIRECT_NZYX_ELEM(img(), z, 0, y, x));
			RawImage<T>::data[(z * this->ydim + y) * this->xdim + x] = value;
		}
	}
}

template <class T> 
RawImage<T> BufferedImage<T>::getRef()
{
	return RawImage<T>(
				RawImage<T>::xdim, 
				RawImage<T>::ydim, 
				RawImage<T>::zdim, 
				RawImage<T>::data);
}


// arithmetical operators

template <class T> inline
BufferedImage<T> operator + (const RawImage<T>& i0, const RawImage<T>& i1)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] + i1.data[i];
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator + (const RawImage<T>& i0, T t)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] + t;
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator - (const RawImage<T>& i0, const RawImage<T>& i1)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] - i1.data[i];
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator - (const RawImage<T>& i0, T t)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] - t;
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator - (const RawImage<T>& i0)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = -i0.data[i];
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator * (const RawImage<T>& i0, const RawImage<T>& i1)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] * i1.data[i];
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator * (const RawImage<T>& i0, T t)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] * t;
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator * (T t, const RawImage<T>& i0)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] * t;
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator / (const RawImage<T>& i0, const RawImage<T>& i1)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] / i1.data[i];
	}
	
	return out;
}

template <class T> inline
BufferedImage<T> operator / (const RawImage<T>& i0, T t)
{
	BufferedImage<T> out(i0.xdim, i0.ydim, i0.zdim);
	
	const size_t s = i0.xdim * i0.ydim * i0.zdim;
	
	for (size_t i = 0; i < s; i++)
	{
		out.data[i] = i0.data[i] / t;
	}
	
	return out;
}


#endif
