/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef TOMO_STACK_H
#define TOMO_STACK_H

#include <string>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/error.h>
#include "projection_IO.h"

template <typename T>
class TomoStack
{
    public:

        TomoStack(){}
		
        TomoStack(std::string imagesFn, int imgCount, std::string angles, 
				  std::string affineTransforms, std::string ctfPath,
                  double angpix, double scaleFactor = 1.0, bool loadImgs = true);
		
		
			std::vector<BufferedImage<T> > images;
			std::vector<gravis::d4Matrix> affineXforms;
			std::vector<gravis::d4Matrix> tiltProjs;
			std::vector<gravis::d4Matrix> worldToImage;
			//std::vector<CTF> ctfs;
			double angpix, scaleFactor;
			
			
		TomoStack<T> extractSubStack(gravis::d3Vector center, int w, int h) const;
			
        TomoStack<T> downsample(int factor, int threads = 1, bool downsampleData = true) const;
		
		void saveImages(std::string path, int f0 = 0, int fc = -1) const;
		void saveLike(std::string fnOut, std::string fnExample) const;

        //std::vector<std::vector<gravis::d2Vector> > loadFiducials(std::string file, double scale);
};


template <typename T>
TomoStack<T> :: TomoStack(std::string imagesFn, int imgCount, std::string angles, 
					   std::string affineTransforms, std::string ctfPath,
                       double angpix, double scaleFactor, bool loadImgs)
:   angpix(angpix),
    scaleFactor(scaleFactor)
{
    size_t ast = imagesFn.find_first_of('*');
    if (ast == std::string::npos)
    {
        BufferedImage<T> stack;
		stack.read(imagesFn);
		
		const int w = stack.xdim;
		const int h = stack.ydim;
		
		if (imgCount < 0) imgCount = stack.zdim;
		
		images.resize(imgCount);
		
		for (size_t i = 0; i < imgCount; i++)
		{
			images[i] = BufferedImage<T>(w,h);
			
			for (size_t y = 0; y < h; y++)
			for (size_t x = 0; x < w; x++)
			{
				images[i](x,y) = stack(x,y,i);
			}
		}
    }
	else
	{
		std::string fnBase = imagesFn.substr(0, ast);
		std::string fnEnd = imagesFn.substr(ast+1);
	
		images.resize(imgCount);
	
		for (int i = 0; i < imgCount; i++)
		{
			std::stringstream sts;
			sts << i;
			std::string fn;
			sts >> fn;
	
			std::string fnn = fnBase+fn+fnEnd;
			std::cout << "reading: " << fnn << "\n";
			images[i].read(fnn);
	
			if (!loadImgs) break;
		}
	}

	gravis::d2Vector center;
    center.x = (images[0].xdim - 1.0)/(2.0 * scaleFactor);
    center.y = images[0].ydim/(2.0 * scaleFactor);

    tiltProjs = ProjectionIO::loadTiltProjections(angles, center.x, center.y);
    affineXforms = ProjectionIO::loadAffineTransforms(affineTransforms, center.x, center.y);
    //ctfs = CtfHelper::loadCtffind4(ctfPath, imgCount, 300.0, 2.7, 0.07);

    if (tiltProjs.size() < imgCount)
    {
        REPORT_ERROR("BackprojectionHelper::backproject: not enough angles in "+angles+".");
    }

    if (affineXforms.size() < imgCount)
    {
        REPORT_ERROR("BackprojectionHelper::backproject: not enough affine transforms in "+affineTransforms+".");
    }

    worldToImage.resize(imgCount);

    for (int i = 0; i < imgCount; i++)
    {
        gravis::d4Matrix Ai = affineXforms[i];
        Ai.invert();

        worldToImage[i] = Ai * tiltProjs[i];

        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 4; k++)
        {
            worldToImage[i](j,k) *= scaleFactor;
        }
    }
}

template <typename T>
TomoStack<T> TomoStack<T> :: extractSubStack(gravis::d3Vector center, int w, int h) const
{
    const int ic = images.size();

    TomoStack<T> ts;

    ts.angpix = angpix;
    ts.scaleFactor = scaleFactor;

    ts.images.resize(ic);
    ts.affineXforms.resize(ic);
    ts.tiltProjs.resize(ic);
    ts.worldToImage.resize(ic);
    //ts.ctfs.resize(ic);

    gravis::d4Vector pw(center.x, center.y, center.z, 1.0);

    for (int i = 0; i < ic; i++)
    {
        gravis::d4Vector pi = worldToImage[i] * pw;
        int x0 = (int)(pi.x - w/2.0 + 0.5);
        int y0 = (int)(pi.y - h/2.0 + 0.5);

        Cutting::extract2D(images[i], ts.images[i], x0, y0, w, h);

        ts.affineXforms[i] = affineXforms[i];
        ts.affineXforms[i](0,3) -= x0;
        ts.affineXforms[i](1,3) -= y0;

        ts.tiltProjs[i] = tiltProjs[i];

        ts.worldToImage[i] = worldToImage[i];
        ts.worldToImage[i](0,3) -= x0;
        ts.worldToImage[i](1,3) -= y0;

        //ts.ctfs[i] = ctfs[i];
    }

    return ts;
}

template <typename T>
TomoStack<T> TomoStack<T> :: downsample(int factor, int threads, bool downsampleData) const
{
	TomoStack<T> out = *this;
	
    const int ic = images.size();
	
	const int w1 = images[0].xdim / factor;
	const int h1 = images[0].ydim / factor;
	
	#pragma omp parallel for num_threads(threads)
    for (int i = 0; i < ic; i++)
    {
		if (downsampleData)
		{
			BufferedImage<T> temp = images[i];
			out.images[i] = Resampling::downsampleFilt_2D_full(temp, w1, h1);
		}
		else
		{
			out.images[i].resize(w1,h1);
		}
		
		out.worldToImage[i] /= factor;
		out.worldToImage[i](3,3) = 1.0;
    }

    out.angpix *= factor;
    out.scaleFactor /= factor;
	
	return out;
}

template <typename T>
void TomoStack<T> :: saveImages(std::string path, int f0, int fc) const
{
    size_t ast = path.find_first_of('*');
	
    if (ast == std::string::npos)
    {
		BufferedImage<T> out(images[0].xdim, images[0].ydim, images.size());
		
		for (size_t n = 0; n < out.zdim; n++)
		for (size_t y = 0; y < out.ydim; y++)
		for (size_t x = 0; x < out.xdim; x++)
		{
			out(x,y,n) = images[n](x,y);
		}
		
		out.write(path);
    }
	else
	{
		std::string fnBase = path.substr(0, ast);
		std::string fnEnd = path.substr(ast+1);
	
		const int ic = fc < 0? images.size() : fc+f0;
	
		for (int i = f0; i < ic; i++)
		{
			std::stringstream sts;
			sts << i;
			std::string fn;
			sts >> fn;
	
			std::string fnn = fnBase+fn+fnEnd;
			std::cout << "writing: " << fnn << "\n";
			images[i].write(fnn);
		}
	}
}

template <typename T>
void TomoStack<T> :: saveLike(std::string fnOut, std::string fnExample) const
{
	const int w = images[0].xdim;
	const int h = images[0].ydim;
	const int fc = images.size();
	
	Image<T> example;
	example.read(fnExample, false);
	
	example.data = MultidimArray<T>(fc, h, w);
	
	for (int f = 0; f < fc; f++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		example.data(f,y,x) = images[f](x,y);
	}
	
	example.write(fnOut);
}

#endif
