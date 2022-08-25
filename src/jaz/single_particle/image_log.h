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

#ifndef IMAGE_LOG_H
#define IMAGE_LOG_H

#include <src/image.h>
#include <src/jaz/single_particle/jaz_config.h>
#include <src/jaz/single_particle/vtk_helper.h>


class ImageLog
{
    public:

		enum Centering {NoCenter, CenterXY, CenterXYZ};


        template <typename T>
        static void write(
            Image<T>& img, std::string fn, Centering center = NoCenter,
            double originX = 0.0, double originY = 0.0, double originZ = 0.0,
            double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0);

        static void write(
            Image<Complex>& img, std::string fn, bool polar, Centering center = NoCenter,
            double originX = 0.0, double originY = 0.0, double originZ = 0.0,
            double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0);

        template <typename T>
        static void write(
            MultidimArray<T>& mda, std::string fn, Centering center = NoCenter,
            double originX = 0.0, double originY = 0.0, double originZ = 0.0,
            double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0);

        template <typename T>
        static void write(
            std::vector<Image<T>>& vec, std::string fn, Centering center = NoCenter,
            double originX = 0.0, double originY = 0.0, double originZ = 0.0,
            double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0);
};

template <typename T>
void ImageLog::write(
        Image<T> &img, std::string fn, Centering center,
        double originX, double originY, double originZ,
        double spacingX, double spacingY, double spacingZ)
{
    const int w = img.data.xdim;
    const int h = img.data.ydim;
    const int d = img.data.zdim;

    switch (center)
    {
        case NoCenter:
	    break;

        case CenterXY:
        {
            Image<T> img2(w,h,d);

            for (int z = 0; z < d; z++)
            for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
            {
                const int xx = (x+w/2)%w;
                const int yy = (y+h/2)%h;

                img2(z,y,x) = img(z,yy,xx);
            }

            write(img2, fn, NoCenter, originX, originY, originZ, spacingX, spacingY, spacingZ);
            return;
        }

        case CenterXYZ:
        {
            Image<T> img2(w,h,d);

            for (int z = 0; z < d; z++)
            for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
            {
                const int xx = (x+w/2)%w;
                const int yy = (y+h/2)%h;
                const int zz = (z+d/2)%d;

                img2(z,y,x) = img(zz,yy,xx);
            }

            write(img2, fn, NoCenter, originX, originY, originZ, spacingX, spacingY, spacingZ);
            return;
        }
    }

    if (JazConfig::writeMrc)
    {
        img.write(fn+".mrc");
    }

    if (JazConfig::writeVtk)
    {
        VtkHelper::writeVTK(img, fn+".vtk", originX, originY, originZ, spacingX, spacingY, spacingZ);
    }
}

template <typename T>
void ImageLog::write(
        MultidimArray<T>& mda, std::string fn, Centering center,
        double originX, double originY, double originZ,
        double spacingX, double spacingY, double spacingZ)
{
    Image<T> img;
    img.data = mda;

    write(img, fn, center, originX, originY, originZ, spacingX, spacingY, spacingZ);
}

template <typename T>
void ImageLog::write(
        std::vector<Image<T>>& vec, std::string fn, Centering center,
        double originX, double originY, double originZ,
        double spacingX, double spacingY, double spacingZ)
{
    if (vec.size() == 0)
    {
        std::cerr << "WARNING: nothing to write to " << fn << " - vector size is zero.\n";
        return;
    }

    if (vec[0].data.zdim > 1 || vec[0].data.ndim > 1)
    {
        REPORT_ERROR("ImageLog::write: unable to write a vector of 3D images\n");
    }

    if (center == CenterXYZ)
    {
        REPORT_ERROR("ImageLog::write: unable to XYZ-center a vector of 2D images\n");
    }

    const int w = vec[0].data.xdim;
    const int h = vec[0].data.ydim;
    const int ic = vec.size();

    Image<T> img(w,h,ic);

    for (int i = 0; i < ic; i++)
    {
        if (vec[i].data.xdim != w || vec[i].data.ydim != h)
        {
            REPORT_ERROR("ImageLog::write: images in vector are of unequal size\n");
        }

        for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            if (center == CenterXY)
            {
                const int xx = (x + w/2) % w;
                const int yy = (y + h/2) % h;

                img(i,y,x) = vec[i](yy,xx);
            }
            else // center == NoCenter
            {
                img(i,y,x) = vec[i](y,x);
            }
        }
    }

    write(img, fn, NoCenter, originX, originY, originZ, spacingX, spacingY, spacingZ);
}

#endif
