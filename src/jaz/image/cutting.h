#ifndef IMG_CUTTING_H
#define IMG_CUTTING_H

#include "buffered_image.h"

class Cutting
{
    public:
		
		template <class T>
		static void extract2D(
				const RawImage<T>& src, 
				BufferedImage<T>& dest,
				long int x0, long int y0,
				long int w, long int h);

		template <class T>
		static BufferedImage<T> extractAxialSlice(
				const RawImage<T>& src,
				int axis,
				int index);
};

template <class T>
void Cutting::extract2D(
		const RawImage<T>& src, 
		BufferedImage<T>& dest,
		long int x0, long int y0,
		long int w, long int h)
{
	const long int d = src.zdim;
	
    dest = BufferedImage<T>(w,h,d);
	
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        long int xx = x0 + x;
        long int yy = y0 + y;

        if (   xx >= 0 && xx < src.xdim
            && yy >= 0 && yy < src.ydim)
        {
			dest(x,y,z) = src(xx,yy,z);
        }
        else
        {
			dest(x,y,z) = T(0);
        }
    }
}

template <class T>
BufferedImage<T> Cutting::extractAxialSlice(
		const RawImage<T>& src,
		int axis,
		int index)
{
	std::vector<long int> dims = src.getSizeVector();

	const int d0 = axis;
	const int d1 = (d0 + 1) % 3;
	const int d2 = (d0 + 2) % 3;

	if (dims[d0] < index)
	{
		REPORT_ERROR_STR("Cutting::extractAxialSlice: bad index: " << index);
	}

	const int w = dims[d1];
	const int h = dims[d2];

	BufferedImage<T> out(w,h);

	std::vector<int> coords(3);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		coords[d0] = index;
		coords[d1] = x;
		coords[d2] = y;

		out(x,y) = src(coords[0], coords[1], coords[2]);
	}

	return out;
}

#endif
