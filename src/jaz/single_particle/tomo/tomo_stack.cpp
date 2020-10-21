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

#include "tomo_stack.h"
#include "projection_helper.h"
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/ctf_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/error.h>

using namespace gravis;

TomoStack :: TomoStack(std::string imagesFn, int imgCount, std::string angles, 
					   std::string affineTransforms, std::string ctfPath,
                       double angpix, double scaleFactor, bool loadImgs)
:   angpix(angpix),
    scaleFactor(scaleFactor)
{
    size_t ast = imagesFn.find_first_of('*');
    if (ast == std::string::npos)
    {
        REPORT_ERROR("TomoStack::ctor: asterisk required in image filename.\n");
    }
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

	d2Vector center;
    center.x = images[0].data.xdim/(2.0 * scaleFactor);
    center.y = images[0].data.ydim/(2.0 * scaleFactor);

    tiltProjs = ProjectionHelper::loadTiltProjections(angles, center.x, center.y);
    affineXforms = ProjectionHelper::loadAffineTransforms(affineTransforms, center.x, center.y);
    ctfs = CtfHelper::loadCtffind4(ctfPath, imgCount, 300.0, 2.7, 0.07);

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
        d4Matrix Ai = affineXforms[i];
        Ai.invert();

        worldToImage[i] = Ai * tiltProjs[i];

        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 4; k++)
        {
            worldToImage[i](j,k) *= scaleFactor;
        }
    }
}

TomoStack TomoStack :: extractSubStack(gravis::d3Vector center, int w, int h)
{
    const int ic = images.size();

    TomoStack ts;

    ts.angpix = angpix;
    ts.scaleFactor = scaleFactor;

    ts.images.resize(ic);
    ts.affineXforms.resize(ic);
    ts.tiltProjs.resize(ic);
    ts.worldToImage.resize(ic);
    ts.ctfs.resize(ic);

    d4Vector pw(center.x, center.y, center.z, 1.0);

    for (int i = 0; i < ic; i++)
    {
        d4Vector pi = worldToImage[i] * pw;
        int x0 = (int)(pi.x - w/2.0 + 0.5);
        int y0 = (int)(pi.y - h/2.0 + 0.5);

        FilterHelper::extract2D(images[i], ts.images[i], x0, y0, w, h);

        ts.affineXforms[i] = affineXforms[i];
        ts.affineXforms[i](0,3) -= x0;
        ts.affineXforms[i](1,3) -= y0;

        ts.tiltProjs[i] = tiltProjs[i];

        ts.worldToImage[i] = worldToImage[i];
        ts.worldToImage[i](0,3) -= x0;
        ts.worldToImage[i](1,3) -= y0;

        ts.ctfs[i] = ctfs[i];
    }

    return ts;
}

void TomoStack :: downsample(int factor, int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    Image<RFLOAT> temp(images[0].data.xdim/factor, images[0].data.ydim/factor);

    for (int i = f0; i < ic; i++)
    {
        SliceHelper::downsample(images[i], temp);
        images[i] = temp;
		
		worldToImage[i] /= factor;
		worldToImage[i](3,3) = 1.0;
    }

    angpix *= factor;
    scaleFactor /= factor;
}

void TomoStack :: phaseFlip(int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    for (int i = f0; i < ic; i++)
    {
        FilterHelper::phaseFlip(images[i], ctfs[i], angpix, images[i]);
    }
}

void TomoStack :: ctfModulate(int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    for (int i = f0; i < ic; i++)
    {
        FilterHelper::modulate(images[i], ctfs[i], angpix, images[i]);
    }
}

void TomoStack :: wienerFilter(RFLOAT eps, RFLOAT Bfac, int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    for (int i = f0; i < ic; i++)
    {
        FilterHelper::wienerFilter(images[i], ctfs[i], angpix, eps, Bfac, images[i]);
    }
}

void TomoStack :: richardsonLucy(int iterations, RFLOAT eps, int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    for (int i = f0; i < ic; i++)
    {
        FilterHelper::richardsonLucy(images[i], ctfs[i], angpix, eps, iterations, images[i]);
    }
}

void TomoStack :: rampFilter(RFLOAT s0, RFLOAT t1, int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    for (int i = f0; i < ic; i++)
    {
        d2Vector u(affineXforms[i](0,0), affineXforms[i](1,0));
        u + u / u.length();

        FilterHelper::rampFilter(images[i], s0, t1, u.x, u.y, images[i]);
    }
}

void TomoStack :: safeLog(RFLOAT eps, int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    for (int im = f0; im < ic; im++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(images[im].data)
        {
            const double v = DIRECT_A2D_ELEM(images[im].data, i, j);

            DIRECT_A2D_ELEM(images[im].data, i, j) = v > eps ? log(v) : log(eps);
        }
    }
}

void TomoStack :: scaledExp(RFLOAT scale, int f0, int fc)
{
    const int ic = fc < 0? images.size() : fc+f0;

    for (int im = f0; im < ic; im++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(images[im].data)
        {
            DIRECT_A2D_ELEM(images[im].data, i, j) = exp(scale * DIRECT_A2D_ELEM(images[im].data, i, j));
        }
    }
}

void TomoStack :: defocusStack(int f, double dz0, double dz1, double eps, double Bfac, std::vector<Image<RFLOAT> >& dest, int x0, int y0, int w, int h)
{
    int zc = dest.size();

    CTF ctf0 = ctfs[f];
    CTF ctf = ctf0;

    Image<RFLOAT> img;

    if (w < 0 || h < 0)
    {
        img = images[f];
    }
    else
    {
        img = Image<RFLOAT>(w,h,1,1);
        FilterHelper::extract2D(images[f], img, x0, y0, w, h);
    }

    for (int z = 0; z < zc; z++)
    {
        double dz = dz0 + z * (dz1 - dz0) / (double)(zc - 1);

        std::cout << z << ": " << dz << " \t ";

        ctf.DeltafU = ctf0.DeltafU + dz;
        ctf.DeltafV = ctf0.DeltafV + dz;

        ctf.initialise();

        FilterHelper::wienerFilter(img, ctf, angpix, eps, Bfac, dest[z]);
        //FilterHelper::lowPassFilter(dest[z], 0.2, 0.1);
        //FilterHelper::phaseFlip(img, ctf, angpix, dest[z]);

        std::cout << FilterHelper::totalLogVariation(dest[z]) << "\n";
    }
}

void TomoStack :: saveImages(std::string path, int f0, int fc)
{
    size_t ast = path.find_first_of('*');
    if (ast == std::string::npos)
    {
        REPORT_ERROR("TomoStack::saveImages: asterisk required in path.\n");
    }

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

std::vector<std::vector<gravis::d2Vector> > TomoStack::loadFiducials(std::string file, double scale)
{
    std::ifstream is(file);

    if (!is.is_open())
    {
        REPORT_ERROR("failed to open " + file + '\n');
    }

    std::vector<d2Vector> batch(0);
    std::vector<std::vector<d2Vector> > out(0);

    int lastF = -1;

    while (is.good())
    {
        double x, y;
        int f;

        is >> x;
        is >> y;
        is >> f;

        if (f < 0 || f >= affineXforms.size())
        {
            std::stringstream sts;
            sts << f;
            std::string fs;
            sts >> fs;

            REPORT_ERROR("illegal fiducial index: "+fs+"\n");
        }

        while (f > lastF + 1)
        {
            std::cout << "Warning: fiducial position for frame " << (lastF + 1) << " is missing.\n";
            lastF++;
            batch.push_back(d2Vector(-1000.0,-1000.0));
        }

        if (f == 0)
        {
            if (batch.size() > 0)
            {
                out.push_back(batch);
            }

            batch.clear();
        }

        lastF = f;

        d4Matrix Ai = affineXforms[f];
        Ai.invert();

        d4Vector d0(scale*x,scale*y,0,1);
        d4Vector d = Ai * d0;

        batch.push_back(d2Vector(d.x,d.y));
    }

    out.push_back(batch);

    return out;
}
