#ifndef NEW_VTK_HELPER_H
#define NEW_VTK_HELPER_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/math/tensor3x3.h>

class NewVtkHelper
{
	public:
		
		template <class T>
		static void writeR3(const RawImage<gravis::t3Vector<T>>& img, std::string fn);
		
		template <class T>
		static void writeTensors(const RawImage<Tensor3x3<T>>& img, std::string fn);
};

template<class T>
void NewVtkHelper::writeR3(
	const RawImage<gravis::t3Vector<T> > &img, std::string fn)
{
	const size_t size = (img.xdim * img.ydim * img.zdim);
    std::ofstream os(fn.c_str());

    os << "# vtk DataFile Version 2.0\n";
    os << "Dynamite(TM) volume\n";
	os << "ASCII\n";
    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.xdim << " " << img.ydim << " " << img.zdim << "\n";
    os << "SPACING 1 1 1\n";
    os << "ORIGIN 0 0 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars double 3\n";
    os << "LOOKUP_TABLE default\n";
	
	for (int z = 0; z < img.zdim; z++)
	for (int y = 0; y < img.ydim; y++)
	for (int x = 0; x < img.xdim; x++)
	{
		os << img(x,y,z).x << " " << img(x,y,z).y << " " << img(x,y,z).z << "\n";
	}
}

template <class T>
void NewVtkHelper::writeTensors(const RawImage<Tensor3x3<T>>& img, std::string fn)
{
	const size_t size = (img.xdim * img.ydim * img.zdim);
    std::ofstream os(fn.c_str());

    os << "# vtk DataFile Version 2.0\n";
    os << "Dynamite(TM) volume\n";
	os << "ASCII\n";
    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.xdim << " " << img.ydim << " " << img.zdim << "\n";
    os << "SPACING 1 1 1\n";
    os << "ORIGIN 0 0 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars double 6\n";
    os << "LOOKUP_TABLE default\n";

	for (int z = 0; z < img.zdim; z++)
	for (int y = 0; y < img.ydim; y++)
	for (int x = 0; x < img.xdim; x++)
	{
		Tensor3x3<T> t = img(x,y,z);
		
		os << t.xx << " " 
		   << t.yy << " " 
		   << t.zz << " " 
		   << t.xy << " " 
		   << t.yz << " " 
		   << t.xz << "\n";
    }
}

#endif
