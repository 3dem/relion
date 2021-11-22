#ifndef JAZ_RAW_IMAGE_H
#define JAZ_RAW_IMAGE_H

#include <vector>
#include <string>
#include <src/macros.h>
#include <src/jaz/gravis/tImage.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/image.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <stddef.h>

#define RAW_IMAGE_BOUNDS_CHECKING 0

template <class T>
class RawImage
{
	public:
		
		static gravis::t3Vector<size_t> readSize(std::string filename);
		
		RawImage()
			: xdim(0), ydim(0), zdim(0)
		{}
		
		RawImage(size_t xdim, size_t ydim, size_t zdim, T* data);
		
		RawImage(Image<T>& img);
		
		RawImage(MultidimArray<T>& mda);
		
			
			long int xdim, ydim, zdim;
			T* data;
		
		
		/* operator (x,y,z): returns a reference to the indicated voxel.
		   The correct version (const or non-const) will be chosen by the compiler,
		   depending on whether the instance is declared const.*/
		
		inline const T& operator() (size_t, size_t, size_t) const;
		inline T& operator() (size_t, size_t, size_t);
		
		inline const T& operator() (size_t, size_t) const;
		inline T& operator() (size_t, size_t);
		
		inline const T& operator[] (size_t) const;
		inline T& operator[] (size_t);
		
		inline const T* getData() const;
		inline T* getData();
		
		inline size_t getSize() const;
		
		std::string getSizeString() const;

		std::vector<long int> getSizeVector() const;
		
		void fill(T t);
		
		RawImage<T> getFullRef();
		RawImage<T> getSliceRef(size_t z);
		const RawImage<T> getConstSliceRef(size_t z) const;
		RawImage<T> getSlabRef(size_t z, size_t thickness);
		const RawImage<T> getConstSlabRef(size_t z, size_t thickness);
		
		template <class T2>
		void copySliceFrom(size_t z, const RawImage<T2>& src, size_t z_src = 0);
		
		void swapWith(RawImage<T>& img);
		
		void write(std::string fn) const;
		void write(std::string fn, double pixelSize, bool writeFloat16 = false) const;
		void writePng(std::string fn) const;
		void writeVtk(
				std::string fn,
				gravis::d3Vector origin = gravis::d3Vector(0,0,0),
				gravis::d3Vector step = gravis::d3Vector(1,1,1)) const;
		void writeNoVTK(std::string fn, double pixelSize) const;
		
		template <class T2>
		void copyTo(gravis::tImage<T2>& img, int z = 0) const;

		template <class T2>
		void copyTo(Image<T2>& img) const;

		template <class T2>
		void copyToN(Image<T2>& img) const;

		template <class T2>
		void copyTo(MultidimArray<T2>& img) const;

		template <class T2>
		void copyToN(MultidimArray<T2>& img) const;
		
		template <class T2>
		void copyTo(RawImage<T2>& img) const;
		
		template <class T2>
		void copyFrom(const Image<T2>& img);
		
		template <class T2>		
		void copyFrom(const RawImage<T2>& img);
		
		bool hasSize(long int w, long int h, long int d);
		bool hasEqualSize(const RawImage<T>& img);
		bool isSquare();
		bool isCubical();
		
		template<class T2>
		bool containsPoint(gravis::t2Vector<T2> point) const;
		
		template<class T2>
		bool containsPoint(gravis::t3Vector<T2> point) const;

		
		template<class T2> inline
		RawImage& operator += (const RawImage<T2>& v)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] += v.data[i];
			}
			
			return *this;
		}
		
		template<class T2> inline
		RawImage& operator -= (const RawImage<T2>& v)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] -= v.data[i];
			}
			
			return *this;
		}
		
		RawImage& operator += (T t)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] += t;
			}
			
			return *this;
		}
		
		RawImage& operator -= (T t)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] -= t;
			}
			
			return *this;
		}
		
		template<class T2> inline
		RawImage& operator *= (const RawImage<T2>& v)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] *= v.data[i];
			}
			
			return *this;
		}
		
		RawImage& operator *= (T t)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] *= t;
			}
			
			return *this;
		}
		
		template<class T2> inline
		RawImage& operator /= (const RawImage<T2>& v)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] /= v.data[i];
			}
			
			return *this;
		}
		
		RawImage& operator /= (T t)
		{
			const long int s = xdim * ydim * zdim;
			
			for (long int i = 0; i < s; i++)
			{
				data[i] /= t;
			}
			
			return *this;
		}

		void addMultiple(T scale, const RawImage<T>& image);

		void addMultiple(const RawImage<T>& scale, const RawImage<T>& image);
};

template <class T>
gravis::t3Vector<size_t> RawImage<T>::readSize(std::string filename)
{
	Image<T> img;
	img.read(filename, false);
	
	return gravis::t3Vector<size_t>(img.data.xdim, img.data.ydim, img.data.zdim);
}

// constructor(s)

template <class T>
RawImage<T>::RawImage(size_t xdim, size_t ydim, size_t zdim, T* data)
	:   xdim(xdim), ydim(ydim), zdim(zdim), data(data)
{
}

template<class T>
RawImage<T>::RawImage(Image<T> &img)
	:   xdim(img.data.xdim), ydim(img.data.ydim), zdim(img.data.zdim), 
		data(img.data.data)
{	
}

template<class T>
RawImage<T>::RawImage(MultidimArray<T> &mda)
	:   xdim(mda.xdim), ydim(mda.ydim), zdim(mda.zdim), 
		data(mda.data)
{	
}

template<class T>
RawImage<T> RawImage<T>::getFullRef()
{
	return RawImage<T>(xdim, ydim, zdim, data);
}

template<class T>
RawImage<T> RawImage<T>::getSliceRef(size_t z)
{
	return RawImage<T>(xdim, ydim, 1, data + xdim * ydim * z);
}

template<class T>
const RawImage<T> RawImage<T>::getConstSliceRef(size_t z) const
{
	return RawImage<T>(xdim, ydim, 1, data + xdim * ydim * z);
}

template<class T>
RawImage<T> RawImage<T>::getSlabRef(size_t z, size_t thickness)
{
	return RawImage<T>(xdim, ydim, thickness, data + xdim * ydim * z);
}

template<class T>
const RawImage<T> RawImage<T>::getConstSlabRef(size_t z, size_t thickness)
{
	return RawImage<T>(xdim, ydim, thickness, data + xdim * ydim * z);
}

template<class T> template<class T2>
void RawImage<T>::copySliceFrom(size_t z_dest, const RawImage<T2>& src, size_t z_src)
{
	for (long int y = 0; y < ydim; y++)
	for (long int x = 0; x < xdim; x++)
	{
		data[z_dest * ydim * xdim + y * xdim + x] = (T) src(x,y,z_src);
	}
}

template<class T>
void RawImage<T>::swapWith(RawImage<T>& img)
{
	T* sw = img.data;
	img.data = data;
	data = sw;
}

// access operators

template <class T> 
inline const T& RawImage<T>::operator() (size_t x, size_t y, size_t z) const
{
	#if RAW_IMAGE_BOUNDS_CHECKING

	const long int i = (z*ydim + y)*xdim + x;

	if (i < 0 || i >= xdim*ydim*zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::operator(): indices "
			<< x << ", " << y << ", " << z
			<< " (" << i << ") not in [0, " << (xdim*ydim*zdim) << "]");
	}

	#endif

	return data[(z*ydim + y)*xdim + x];
}

template <class T>
inline T& RawImage<T>::operator() (size_t x, size_t y, size_t z)
{
	#if RAW_IMAGE_BOUNDS_CHECKING

	const long int i = (z*ydim + y)*xdim + x;

	if (i < 0 || i >= xdim*ydim*zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::operator(): indices "
			<< x << ", " << y << ", " << z
			<< " (" << i << ") not in [0, " << (xdim*ydim*zdim) << "]");
	}

	#endif

	return data[(z*ydim + y)*xdim + x];
}

template <class T>
inline const T& RawImage<T>::operator() (size_t x, size_t y) const
{
	#if RAW_IMAGE_BOUNDS_CHECKING

	const long int i = y*xdim + x;

	if (i < 0 || i >= xdim*ydim*zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::operator(): indices "
			<< x << ", " << y
			<< " (" << i << ") not in [0, " << (xdim*ydim*zdim) << "]");
	}

	#endif

	return data[y*xdim + x];
}

template <class T>
inline T& RawImage<T>::operator() (size_t x, size_t y)
{
	#if RAW_IMAGE_BOUNDS_CHECKING

	const long int i = y*xdim + x;

	if (i < 0 || i >= xdim*ydim*zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::operator(): indices "
			<< x << ", " << y
			<< " (" << i << ") not in [0, " << (xdim*ydim*zdim) << "]");
	}

	#endif


	return data[y*xdim + x];
}

template <class T>
inline const T& RawImage<T>::operator[] (size_t i) const
{
	#if RAW_IMAGE_BOUNDS_CHECKING

	if (i >= xdim*ydim*zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::operator(): index "
			<< i << " not in [0, " << (xdim*ydim*zdim) << "]");
	}

	#endif

	return data[i];
}

template <class T>
inline T& RawImage<T>::operator[] (size_t i)
{
	#if RAW_IMAGE_BOUNDS_CHECKING

	if (i >= xdim*ydim*zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::operator(): index "
			<< i << " not in [0, " << (xdim*ydim*zdim) << "]");
	}

	#endif

	return data[i];
}

// utility functions

template <class T>
inline const T* RawImage<T>::getData() const
{
	return data;
}

template <class T>
inline T* RawImage<T>::getData()
{
	return data;
}

template <class T>
inline size_t RawImage<T>::getSize() const
{
	return xdim * ydim * zdim;
}

template <class T>
inline std::string RawImage<T>::getSizeString() const
{
	std::ostringstream sts;
	sts << xdim << 'x' << ydim << 'x' << zdim;
	return sts.str();
}

template <class T>
std::vector<long int> RawImage<T>::getSizeVector() const
{
	std::vector<long int> out(3);
	out[0] = xdim;
	out[1] = ydim;
	out[2] = zdim;
	return out;
}

template <class T>
inline void RawImage<T>::fill(T t)
{
	const long int s = xdim * ydim * zdim;
	
	for (long int i = 0; i < s; i++)
	{
		data[i] = t;
	}
}

// I/O

template <class T>
void RawImage<T>::write(std::string fn) const
{
	write(fn, 1);
}

template <class T>
void RawImage<T>::write(std::string fn, double pixelSize, bool writeFloat16) const
{
	std::string::size_type dot = fn.find_last_of('.');
	
	if (dot == std::string::npos)
	{
		REPORT_ERROR_STR("RawImage<T>::write: filename has no ending (" << fn << ").");
	}
	
	std::string end = fn.substr(dot+1);
	
	if (end == "vtk")
	{
		if (writeFloat16)
		{
			REPORT_ERROR_STR("RawImage<T>::write: Float16 type not implemented for VTK format.");
		}
		Image<T> img;
		copyTo(img);
		VtkHelper::writeVTK(img, fn, 0, 0, 0, pixelSize, pixelSize, pixelSize);
	}
	else if (end == "mrcs")
	{
		Image<T> img;
		copyToN(img);
		img.setSamplingRateInHeader(pixelSize);
		img.write(fn);
		img.write(fn, -1, false, WRITE_OVERWRITE, writeFloat16 ? Float16: Float);

	}
	else
	{
		Image<T> img;
		copyTo(img);
		img.setSamplingRateInHeader(pixelSize);
		img.write(fn, -1, false, WRITE_OVERWRITE, writeFloat16 ? Float16: Float);
	}
}


template <class T>
void RawImage<T>::writePng(std::string fn) const
{
	if (zdim > 1) 
	{
		REPORT_ERROR("RawImage<T>::write: cannot write 3D image into PNG.");
	}
	
	gravis::tImage<T> img;
	copyTo(img);
	img.write(fn);
}

template <class T>
void RawImage<T>::writeVtk(std::string fn, gravis::d3Vector origin, gravis::d3Vector step) const
{
	Image<T> img;
	copyTo(img);
	
	VtkHelper::writeVTK(
				img, fn, 
				origin.x, origin.y, origin.z, 
				step.x, step.y, step.z);
}

template<class T>
void RawImage<T>::writeNoVTK(std::string fn, double pixelSize) const
{
	std::string::size_type dot = fn.find_last_of('.');

	if (dot == std::string::npos)
	{
		REPORT_ERROR_STR("RawImage<T>::write: filename has no ending (" << fn << ").");
	}

	std::string end = fn.substr(dot+1);

	if (end == "mrcs")
	{
		Image<T> img;
		copyToN(img);
		img.setSamplingRateInHeader(pixelSize);
		img.write(fn);
	}
	else
	{
		Image<T> img;
		copyTo(img);
		img.setSamplingRateInHeader(pixelSize);
		img.write(fn);
	}
}

template <class T> template <class T2>
inline void RawImage<T>::copyTo(gravis::tImage<T2>& img, int z) const
{
	img = gravis::tImage<T2>(xdim, ydim);
	
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		img(x,y) = (T2) data[(z*ydim + y)*xdim + x];
	}
}

template <class T> template <class T2>
inline void RawImage<T>::copyTo(Image<T2>& img) const
{
	copyTo(img.data);
}

template <class T> template <class T2>
inline void RawImage<T>::copyToN(Image<T2>& img) const
{
	copyToN(img.data);
}

template <class T> template <class T2>
inline void RawImage<T>::copyTo(MultidimArray<T2>& img) const
{
	img = MultidimArray<T2>(zdim, ydim, xdim);

	for (size_t z = 0; z < zdim; z++)
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		DIRECT_A3D_ELEM(img, z, y, x) = (T2) data[(z*ydim + y)*xdim + x];
	}
}

template <class T> template <class T2>
inline void RawImage<T>::copyToN(MultidimArray<T2>& img) const
{
	img = MultidimArray<T2>(zdim, 1, ydim, xdim);

	for (size_t z = 0; z < zdim; z++)
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		DIRECT_NZYX_ELEM(img, z, 0, y, x) = (T2) data[(z*ydim + y)*xdim + x];
	}
}

template <class T> template <class T2>
inline void RawImage<T>::copyTo(RawImage<T2>& img) const
{
	for (size_t z = 0; z < zdim; z++)
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		img(x,y,z) = (T2) data[(z*ydim + y)*xdim + x];
	}
}

template <class T> template <class T2>
inline void RawImage<T>::copyFrom(const Image<T2>& img)
{
	if (xdim > img.data.xdim || ydim > img.data.ydim || zdim > img.data.zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::copyFrom: input image too small: "
			<< img.data.xdim << "x" << img.data.ydim << "x" << img.data.zdim 
			<< " (should be at least: " 
			<< xdim << "x" << ydim << "x" << zdim << ")");
	}
	
	for (size_t z = 0; z < zdim; z++)
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		T val = (T) DIRECT_A3D_ELEM(img(), z, y, x);
		data[(z*ydim + y)*xdim + x] = val;
	}
}

template <class T> template <class T2>
inline void RawImage<T>::copyFrom(const RawImage<T2>& img)
{
	if (xdim > img.xdim || ydim > img.ydim || zdim > img.zdim)
	{
		REPORT_ERROR_STR("RawImage<T>::copyFrom: input image too small: "
			<< img.xdim << "x" << img.ydim << "x" << img.zdim 
			<< " (should be at least: " 
			<< xdim << "x" << ydim << "x" << zdim << ")");
	}
	
	for (size_t z = 0; z < zdim; z++)
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		data[(z*ydim + y)*xdim + x] = T(img(x,y,z));
	}
}

template <class T>
bool RawImage<T>::hasSize(long int w, long int h, long int d)
{
	return xdim == w && ydim == h && zdim == d;
}

template <class T>
bool RawImage<T>::hasEqualSize(const RawImage<T>& img)
{
	return xdim == img.xdim && ydim == img.ydim && zdim == img.zdim;
}

template<class T>
bool RawImage<T>::isSquare()
{
	return xdim == ydim;
}

template<class T>
bool RawImage<T>::isCubical()
{
	return xdim == ydim && ydim == zdim;
}

template<class T> template<class T2>
bool RawImage<T>::containsPoint(gravis::t2Vector<T2> point) const
{
	return 
			point.x > 0 && point.x < xdim &&
			point.y > 0 && point.y < ydim;
}

template<class T> template<class T2>
bool RawImage<T>::containsPoint(gravis::t3Vector<T2> point) const
{
	return 
			point.x > 0 && point.x < xdim &&
			point.y > 0 && point.y < ydim &&
			point.z > 0 && point.z < zdim;
}

template <class T1, class T2>
RawImage<T1> operator * (const RawImage<T1>& t1, const RawImage<T2>& t2)
{
	RawImage<T1> out = t1;
	out += t2;
	return out;
}

template <class T>
void RawImage<T>::addMultiple(T scale, const RawImage<T>& image)
{
	for (size_t z = 0; z < zdim; z++)
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		data[(z*ydim + y)*xdim + x] += T(scale * image(x,y,z));
	}
}

template <class T>
void RawImage<T>::addMultiple(const RawImage<T>& scale, const RawImage<T>& image)
{
	for (size_t z = 0; z < zdim; z++)
	for (size_t y = 0; y < ydim; y++)
	for (size_t x = 0; x < xdim; x++)
	{
		data[(z*ydim + y)*xdim + x] += T(scale(x,y,z) * image(x,y,z));
	}
}

#endif
