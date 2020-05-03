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

#endif
