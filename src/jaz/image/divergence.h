#ifndef IMAGE_DIVERGENCE_H
#define IMAGE_DIVERGENCE_H

#include "raw_image.h"
#include <src/jaz/gravis/t3Vector.h>

class Divergence
{
	public:
		
		template<typename T>
		static inline void backward3D_inSitu(
				const RawImage<gravis::t3Vector<T>>& img, RawImage<T>& dest);
		
};

template<typename T> 
inline void Divergence::backward3D_inSitu(
		const RawImage<gravis::t3Vector<T>>& img, RawImage<T>& dest)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		T div = 0;
		
		if (x > 0) div += img(x,y,z).x - img(x-1,y,z).x;
		else div += 2 * img(x,y,z).x;
		
		if (y > 0) div += img(x,y,z).y - img(x,y-1,z).y;
		else div += 2 * img(x,y,z).y;
		
		if (z > 0) div += img(x,y,z).z - img(x,y,z-1).z;
		else div += 2 * img(x,y,z).z;
		
		dest(x,y,z) = div;
	}
}

#endif
