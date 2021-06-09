#ifndef IMG_PADDING_H
#define IMG_PADDING_H

#include "raw_image.h"
#include "buffered_image.h"

class Padding
{
	public:
		
		template <class T>
		static BufferedImage<T> padCorner2D_half(const RawImage<T> &img, int w, int h);
		
		template <class T>
		static BufferedImage<T> padCorner2D_full(const RawImage<T> &img, int w, int h);
		
		template <class T>
		static BufferedImage<T> padCenter2D_full(const RawImage<T> &img, int b);
		
		template <class T>
		static BufferedImage<T> padCenter3D_full(const RawImage<T> &img, int b);
		
		template <class T>
		static BufferedImage<T> unpadCenter2D_full(const RawImage<T> &img, int b);
		
		template <class T>
		static BufferedImage<T> unpadCorner2D_full(const RawImage<T> &img, int b);
		
		template <class T>
		static BufferedImage<T> unpadCenter3D_full(const RawImage<T> &img, int b);
		
		template <class T>
		static void copyUnpaddedCenter2D_full(const RawImage<T>& src, RawImage<T>& dest);
};

template <class T>
BufferedImage<T> Padding::padCorner2D_half(const RawImage<T>& img, int w, int h)
{
	const int w0 = img.xdim;
	const int h0 = img.ydim;
	
	BufferedImage<T> out(w,h);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const int x1 = x;
		const int y1 = y < h/2? y : y - h;
		
		if (x1 < w0/2 && y1 < h0/2 && x1 >= -w0/2 && y1 >= -h0/2)
		{
			const int x0 = x1;
			const int y0 = y1 < 0? y1 + h0 : y1;
			
			out(x,y) = img(x0,y0);
		}
		else
		{
			out(x,y) = T(0);
		}
	}
	
	return out;
}

template <class T>
BufferedImage<T> Padding::padCorner2D_full(const RawImage<T>& img, int w, int h)
{
	const int w0 = img.xdim;
	const int h0 = img.ydim;
	
	BufferedImage<T> out(w,h);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const int x1 = x < w/2? x : x - w;
		const int y1 = y < h/2? y : y - h;
		
		if (x1 < w0/2 && y1 < h0/2 && x1 >= -w0/2 && y1 >= -h0/2)
		{
			const int x0 = x1 < 0? x1 + w0 : x1;
			const int y0 = y1 < 0? y1 + h0 : y1;
			
			out(x,y) = img(x0,y0);
		}
		else
		{
			out(x,y) = T(0);
		}
	}
	
	return out;
}

template <class T>
BufferedImage<T> Padding::padCenter2D_full(const RawImage<T> &img, int b)
{
	const int w0 = img.xdim;
	const int h0 = img.ydim;
	
	BufferedImage<T> out(w0 + 2*b, h0 + 2*b, 1);
	out.fill(0);
	
	for (int y = 0; y < h0; y++)
	for (int x = 0; x < w0; x++)
	{
		out(x+b,y+b) = img(x,y);
	}
	
	return out;
}

template <class T>
BufferedImage<T> Padding::padCenter3D_full(const RawImage<T> &img, int b)
{
	const int w0 = img.xdim;
	const int h0 = img.ydim;
	const int d0 = img.zdim;
	
	BufferedImage<T> out(w0 + 2*b, h0 + 2*b, d0 + 2*b);
	out.fill(0);
	
	for (int z = 0; z < d0; z++)
	for (int y = 0; y < h0; y++)
	for (int x = 0; x < w0; x++)
	{
		out(x+b,y+b,z+b) = img(x,y,z);
	}
	
	return out;
}

template <class T>
BufferedImage<T> Padding::unpadCenter2D_full(const RawImage<T> &img, int b)
{
	const int w0 = img.xdim - 2*b;
	const int h0 = img.ydim - 2*b;
	const int fc = img.zdim;
	
	BufferedImage<T> out(w0,h0,fc);
	
	for (int f = 0; f < fc; f++)
	for (int y = 0; y < h0; y++)
	for (int x = 0; x < w0; x++)
	{
		out(x,y,f) = img(x+b, y+b, f);
	}

	return out;
}

template <class T>
BufferedImage<T> Padding::unpadCorner2D_full(const RawImage<T> &img, int b)
{
	const int w0 = img.xdim - 2*b;
	const int h0 = img.ydim - 2*b;
	const int fc = img.zdim;
	
	BufferedImage<T> out(w0,h0,fc);
	
	for (int f = 0; f < fc; f++)
	for (int y = 0; y < h0; y++)
	for (int x = 0; x < w0; x++)
	{
		const int xx = x < w0/2? x : x - w0 + img.xdim;
		const int yy = y < h0/2? y : y - h0 + img.ydim;
		
		out(x,y,f) = img(xx, yy, f);
	}
	
	return out;
}

template <class T>
BufferedImage<T> Padding::unpadCenter3D_full(const RawImage<T> &img, int b)
{
	const int w0 = img.xdim - 2*b;
	const int h0 = img.ydim - 2*b;
	const int d0 = img.zdim - 2*b;
	
	BufferedImage<T> out(w0,h0,d0);
	
	for (int z = 0; z < d0; z++)
	for (int y = 0; y < h0; y++)
	for (int x = 0; x < w0; x++)
	{
		out(x,y,z) = img(x+b, y+b, z+b);
	}
	
	return out;
}

template<class T>
void Padding::copyUnpaddedCenter2D_full(
		const RawImage<T>& src, 
		RawImage<T>& dest)
{
	const int w_src = src.xdim;
	const int h_src = src.ydim;
	
	const int w_dest = dest.xdim;
	const int h_dest = dest.ydim;
	
	const int bx = (w_src - w_dest) / 2;
	const int by = (h_src - h_dest) / 2;
	
	for (int y = 0; y < h_dest; y++)
	for (int x = 0; x < w_dest; x++)
	{
		dest(x,y) = src(x + bx, y + by);
	}
}

#endif
