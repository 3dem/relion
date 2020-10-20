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


#include <src/jaz/gravis/t4Matrix.h>
#include "backprojection_helper.h"
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/config.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/interpolation.h>
#include <src/jaz/single_particle/image_log.h>

#include <src/fftw.h>

using namespace gravis;

void BackprojectionHelper::backprojectRaw(
        const Image<RFLOAT>& stack, std::string tiltAngles, Volume<RFLOAT>& dest, Volume<unsigned char>& maskDest,
        d3Vector origin, double spacing, int frames)
{
    const double cix = stack.data.xdim / 2.0;
    const double ciy = stack.data.ydim / 2.0;

    std::cout << "center: " << cix << ", " << ciy << "\n";

    std::ifstream anglesFile(tiltAngles.c_str());

    if (!anglesFile.is_open())
    {
        REPORT_ERROR("BackprojectionHelper::backproject: failed to open "+tiltAngles+".");
    }

    // vol2world * (0,0,0,1)       = (xw0, yw0, zw0, 1)
    // vol2world * (xwc,ywc,zwc,1) = (xw1, yw1, zw1, 1)

    d4Matrix vol2world;

    vol2world(0,0) = spacing;
    vol2world(1,1) = spacing;
    vol2world(2,2) = spacing;
    vol2world(0,3) = origin.x;
    vol2world(1,3) = origin.y;
    vol2world(2,3) = origin.z;

    std::vector<double> angles;
    std::vector<d4Matrix> vol2img;

    const double deg2rad = PI/180.0;

    while (anglesFile.good())
    {
        double a;
        anglesFile >> a;
        a *= deg2rad;

        angles.push_back(a);

        d4Matrix w2i;
        w2i(0,0) =  cos(a);
        w2i(0,2) =  sin(a);
        w2i(2,0) = -sin(a);
        w2i(2,2) =  cos(a);
        w2i(0,3) =  cix;
        w2i(1,3) =  ciy;

        vol2img.push_back(w2i * vol2world);
    }

    const int ic = frames > 0? frames : stack.data.ndim;

    if (vol2img.size() < ic)
    {
        REPORT_ERROR("BackprojectionHelper::backproject: not enough angles in "+tiltAngles+".");
    }

    #if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif
    FOR_ALL_VOXELS(dest)
    {
        if (x == 0 && y == 0)
        {
            std::cout << z << "/" << dest.dimz << "\n";
        }
        double sum = 0.0;
        double wgh = 0.0;

        d4Vector pw(x,y,z,1.0);

        for (int im = 0; im < ic; im++)
        {
            d4Vector pi = vol2img[im] * pw;

            if (Interpolation::isInSlice(stack, pi.x, pi.y))
            {
                sum += Interpolation::linearXY(stack, pi.x, pi.y, im);
                wgh += 1.0;
            }
        }

        if (wgh > 0.0)
        {
            sum /= wgh;
        }

        dest(x,y,z) = sum;
        maskDest(x,y,z) = wgh > 0.0? 1 : 0;
    }
}

void BackprojectionHelper::backprojectRaw(
		const TomoStack& stack,
        Volume<RFLOAT>& dest, Volume<RFLOAT>& maskDest,
        gravis::d3Vector origin, double spacing, 
		InterpolationType interpolation, double taperX, double taperY, 
		double wMin, int frame0, int frames)
{
    d4Matrix vol2world;

    vol2world(0,0) = spacing;
    vol2world(1,1) = spacing;
    vol2world(2,2) = spacing;
    vol2world(0,3) = origin.x;
    vol2world(1,3) = origin.y;
    vol2world(2,3) = origin.z;
	
	/*std::cout << "vol2world: \n" << vol2world << "\n";
	std::cout << "stack.worldToImage[0]: \n" << stack.worldToImage[0] << "\n";
	std::cout << "vol2img[0]: \n" << (stack.worldToImage[0] * vol2world) << "\n";*/
	

    const int ic = frames > 0? frames + frame0 : stack.images.size();

    std::cout << frame0 << " - " << (ic-1) << "\n";

    std::vector<d4Matrix> vol2img(ic);

    for (int im = 0; im < ic; im++)
    {
        vol2img[im] = stack.worldToImage[im] * vol2world;
    }

    /*#if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif*/
	
	/*std::vector<Image<RFLOAT>> debugImgs(ic);
	for (int im = frame0; im < ic; im++)
	{
		debugImgs[im] = stack.images[im];
	}*/
	
    FOR_ALL_VOXELS(dest)
    {
        double sum = 0.0;
        double wgh = 0.0;

        d4Vector pw(x,y,z,1.0);

        for (int im = frame0; im < ic; im++)
        {
            d4Vector pi = vol2img[im] * pw;

            if (Interpolation::isInSlice(stack.images[im], pi.x, pi.y))
            {
                double wghi = Interpolation::getTaperWeight(stack.images[im], pi.x, pi.y, taperX, taperY);

                if (interpolation == Linear)
                {
                    sum += wghi * Interpolation::linearXY(stack.images[im], pi.x, pi.y, 0);
                }
                else
                {
                    sum += wghi * Interpolation::cubicXY(stack.images[im], pi.x, pi.y, 0);
                }
				
				//debugImgs[im]((int)(pi.y+0.5), (int)(pi.x+0.5)) += 1000.0;

                wgh += wghi;
            }
        }

        if (wgh > 0.0)
        {
            sum /= wgh;
        }

        dest(x,y,z) = sum;
        maskDest(x,y,z) = wgh;
    }
	
	/*JazConfig::writeMrc = false;
	JazConfig::writeVtk = true;
	ImageLog::write(debugImgs, "debug_imgs");*/

    double mean = 0.0, sum = 0.0;

    FOR_ALL_VOXELS(dest)
    {
        mean += maskDest(x,y,z)*dest(x,y,z);
        sum += maskDest(x,y,z);
    }

	if (sum > 0.0)
	{
		mean /= sum;
	}

    #if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif
    FOR_ALL_VOXELS(dest)
    {
        double t = maskDest(x,y,z) / wMin;

        if (t < 1.0)
        {
            dest(x,y,z) = t * dest(x,y,z) + (1.0 - t) * mean;
        }
    }
}


void BackprojectionHelper::backprojectExactWeights(
        const TomoStack& stack, Volume<RFLOAT>& dest,
        d3Vector origin, double spacing, double taperX, double taperY, double taperZ, double wMin,
        int frame0, int frames)
{
    const int wv = dest.dimx;
    const int hv = dest.dimy;
    const int dv = dest.dimz;

    Volume<RFLOAT> vol(wv,hv,dv), volM(wv,hv,dv);

    std::cout << "performing unweighted backprojection...\n";
    backprojectRaw(stack, vol, volM, origin, spacing, Linear, taperX, taperY, wMin, frame0, frames);

    taperEdges(vol, taperX, taperY, taperZ);

    Volume<RFLOAT> weight(wv/2 + 1,hv,dv);

    std::cout << "backprojecting dots...\n";
    backprojectDots(stack, weight, origin, spacing, taperX, taperY, taperZ, frame0, frames);

    std::cout << "applying weights...\n";
    Image<RFLOAT> volRL;
    VolumeConverter::convert(vol, volRL);
    vol.resize(0,0,0);

    Image<Complex> dataFreq;

    FourierTransformer ft;
    ft.FourierTransform(volRL(), dataFreq.data, false);

    FilterHelper::divideExcessive(dataFreq, weight, weight(0,0,0)/(double)stack.images.size(), dataFreq);

    FourierTransformer ft2;
    ft2.inverseFourierTransform(dataFreq.data, volRL());

    VolumeConverter::convert(volRL, dest);
}




void BackprojectionHelper::backprojectExactWeightsFreq(
        const TomoStack& stack,
        Image<Complex>& dest, Volume<RFLOAT>& weight,
        d3Vector origin, double spacing, double taperX, double taperY, double taperZ, double wMin,
        int frame0, int frames)
{
    const int wv = 2 * (dest.data.xdim - 1);
    const int hv = dest.data.ydim;
    const int dv = dest.data.zdim;

    std::cout << wv << "x" << hv << "x" << dv << "x" << dest.data.ndim << "\n";

    Volume<RFLOAT> vol(wv,hv,dv), volM(wv,hv,dv);

    std::cout << "performing unweighted backprojection...\n";
    backprojectRaw(stack, vol, volM, origin, spacing, Linear, taperX, taperY, wMin, frame0, frames);


    taperEdges(vol, taperX, taperY, taperZ);

    weight.resize(wv/2 + 1,hv,dv);

    std::cout << "backprojecting dots...\n";
    backprojectDots(stack, weight, origin, spacing, taperX, taperY, taperZ, frame0, frames);

    std::cout << "applying weights...\n";
    Image<RFLOAT> volRL;
    VolumeConverter::convert(vol, volRL);
    vol.resize(0,0,0);

    FourierTransformer ft;
    ft.FourierTransform(volRL(), dest.data, true);

    const double theta = weight(0,0,0)/(double)stack.images.size();

    for (long int z = 0; z < dv; z++)
    for (long int y = 0; y < hv; y++)
    for (long int x = 0; x < wv/2 + 1; x++)
    {
        const double t = weight(x,y,z)/theta;

        if (t > 1)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) / t;
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(dest.data, 0, z, y, x);
        }

        weight(x,y,z) = t;
    }
}

void BackprojectionHelper::backprojectDots(
    const TomoStack& stack, Volume<RFLOAT>& dest, gravis::d3Vector origin, double spacing, double taperX, double taperY, double taperZ, int frame0, int frames)
{
    const int wv = 2*(dest.dimx-1);
    const int hv = dest.dimy;
    const int dv = dest.dimz;

    d4Matrix vol2world;

    vol2world(0,0) = spacing;
    vol2world(1,1) = spacing;
    vol2world(2,2) = spacing;
    vol2world(0,3) = origin.x;
    vol2world(1,3) = origin.y;
    vol2world(2,3) = origin.z;

    const int ic = frames > 0? frames + frame0 : stack.images.size();

    std::cout << frame0 << " - " << (ic-1) << "\n";

    d4Vector originVol(wv/2, hv/2, dv/2, 1.0);
    std::cout << "originVol = " << originVol << "\n";

    Volume<RFLOAT> streakVol(wv, hv, dv);
    Image<RFLOAT> volRL;
    Image<Complex> spectrum;

    std::vector<d4Matrix> vol2img(ic);
    std::vector<d4Vector> volOrigImg(ic);

    for (int im = 0; im < ic; im++)
    {
        //std::cout << "   " << im << "/" << (ic-1) << "\n";

        vol2img[im] = stack.worldToImage[im] * vol2world;
        volOrigImg[im] = vol2img[im] * originVol;
    }

    dest.fill(0.0);
    streakVol.fill(0.0);

    #if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif
    FOR_ALL_VOXELS(streakVol)
    {
        for (int im = 0; im < ic; im++)
        {
            d4Vector pw(x,y,z,1.0);
            d4Vector pi = vol2img[im] * pw;
            d4Vector d = pi - volOrigImg[im];

            double dx, dy;

            if (d.x < -1 || d.x > 1) dx = 0.0;
            else dx = 1.0 - std::abs(d.x);
            if (d.y < -1 || d.y > 1) dy = 0.0;
            else dy = 1.0 - std::abs(d.y);

            streakVol(x,y,z) += dx*dy;
        }
    }

    taperEdges(streakVol, taperX, taperY, taperZ);
    VolumeConverter::convert(streakVol, volRL);

    CenterFFT(volRL.data, true);

    FourierTransformer ft;
    ft.FourierTransform(volRL(), spectrum.data, false);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(spectrum.data)
    {
        Complex z = DIRECT_A3D_ELEM(spectrum.data, k, i, j);

        dest(j,i,k) = z.abs();
    }
}


void BackprojectionHelper::backprojectDotsFS(
    const TomoStack& stack, Image<Complex>& dest, gravis::d3Vector origin, double spacing, double taperX, double taperY, double taperZ, int frame0, int frames)
{
    const int wv = 2*(dest.data.xdim-1);
    const int hv = dest.data.ydim;
    const int dv = dest.data.zdim;

    d4Matrix vol2world;

    vol2world(0,0) = spacing;
    vol2world(1,1) = spacing;
    vol2world(2,2) = spacing;
    vol2world(0,3) = origin.x;
    vol2world(1,3) = origin.y;
    vol2world(2,3) = origin.z;

    const int ic = frames > 0? frames + frame0 : stack.images.size();

    std::cout << frame0 << " - " << (ic-1) << "\n";

    d4Vector originVol(wv/2, hv/2, dv/2, 1.0);
    std::cout << "originVol = " << originVol << "\n";

    Volume<RFLOAT> streakVol(wv, hv, dv);
    Image<RFLOAT> volRL;

    std::vector<d4Matrix> vol2img(ic);
    std::vector<d4Vector> volOrigImg(ic);

    for (int im = 0; im < ic; im++)
    {
        std::cout << "   " << im << "/" << (ic-1) << "\n";

        vol2img[im] = stack.worldToImage[im] * vol2world;
        volOrigImg[im] = vol2img[im] * originVol;
    }

    streakVol.fill(0.0);

    #if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif
    FOR_ALL_VOXELS(streakVol)
    {
        for (int im = 0; im < ic; im++)
        {
            d4Vector pw(x,y,z,1.0);
            d4Vector pi = vol2img[im] * pw;
            d4Vector d = pi - volOrigImg[im];

            double dx, dy;

            if (d.x < -1 || d.x > 1) dx = 0.0;
            else dx = 1.0 - std::abs(d.x);
            if (d.y < -1 || d.y > 1) dy = 0.0;
            else dy = 1.0 - std::abs(d.y);

            streakVol(x,y,z) += dx*dy;
        }
    }

    taperEdges(streakVol, taperX, taperY, taperZ);
    VolumeConverter::convert(streakVol, volRL);

    CenterFFT(volRL.data, true);

    FourierTransformer ft;
    ft.FourierTransform(volRL(), dest.data, true);
}

void BackprojectionHelper::backprojectDotsSeparately(
    const TomoStack& stack, Volume<RFLOAT>& dest, gravis::d3Vector origin, double spacing, double taperX, double taperY, double taperZ, int frame0, int frames)
{
    const int wv = 2*(dest.dimx-1);
    const int hv = dest.dimy;
    const int dv = dest.dimz;

    d4Matrix vol2world;

    vol2world(0,0) = spacing;
    vol2world(1,1) = spacing;
    vol2world(2,2) = spacing;
    vol2world(0,3) = origin.x;
    vol2world(1,3) = origin.y;
    vol2world(2,3) = origin.z;

    const int ic = frames > 0? frames + frame0 : stack.images.size();

    std::cout << frame0 << " - " << (ic-1) << "\n";

    d4Vector originVol((double)(dest.dimx - 1), (double)dest.dimy/2.0, (double)dest.dimz/2.0, 1.0);

    Volume<RFLOAT> streakVol(wv, hv, dv);
    Image<RFLOAT> volRL;
    Image<Complex> spectrum;

    dest.fill(0.0);

    for (int im = 0; im < ic; im++)
    {
        std::cout << "   " << im << "/" << (ic-1) << "\n";

        d4Matrix vol2img = stack.worldToImage[im] * vol2world;
        d4Vector volOrigImg = vol2img * originVol;

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        FOR_ALL_VOXELS(streakVol)
        {
            double sum = 0.0;

            d4Vector pw(x,y,z,1.0);

            d4Vector pi = vol2img * pw;

            d4Vector d = pi - volOrigImg;

            double dx, dy;

            if (d.x < -1 || d.x > 1) dx = 0.0;
            else dx = 1.0 - std::abs(d.x);
            if (d.y < -1 || d.y > 1) dy = 0.0;
            else dy = 1.0 - std::abs(d.y);

            sum += dx*dy;

            streakVol(x,y,z) = sum;
        }

        taperEdges(streakVol, 0.05*wv, 0.02*hv, 0.05*dv);
        VolumeConverter::convert(streakVol, volRL);

        FourierTransformer ft;
        ft.FourierTransform(volRL(), spectrum.data, false);

        if (im % 10 == 1)
        {
            std::stringstream sts;
            sts << im;
            std::string fn;
            sts >> fn;

            //VtkHelper::writeVTK(streakVol, "streakVol_"+fn+".vtk");
            //VtkHelper::writeVTK_Complex(spectrum.data, "streakVol_"+fn+"_FS.vtk");
        }

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(spectrum.data)
        {
            Complex z = DIRECT_A3D_ELEM(spectrum.data, k, i, j);

            dest(j,i,k) += z.abs();
        }
    }
}



void BackprojectionHelper::backprojectOriginDot(
        const TomoStack& stack,
        Volume<RFLOAT>& dest, double sigma,
        gravis::d3Vector origin, double spacing,
        int frame0, int frames)
{
    d4Matrix vol2world;

    vol2world(0,0) = spacing;
    vol2world(1,1) = spacing;
    vol2world(2,2) = spacing;
    vol2world(0,3) = origin.x;
    vol2world(1,3) = origin.y;
    vol2world(2,3) = origin.z;

    const int ic = frames > 0? frames + frame0 : stack.images.size();

    std::cout << frame0 << " - " << (ic-1) << "\n";

    d4Vector originVol0(0.5/(double)dest.dimx, 0.5/(double)dest.dimy, 0.5/(double)dest.dimz, 1.0);

    std::vector<d4Vector> originVol(8);

    for (int c = 0; c < 8; c++)
    {
        int sx = c%2;
        int sy = (c/2)%2;
        int sz = (c/4)%2;

        originVol[c].x = originVol0.x + sx * (double)dest.dimx;
        originVol[c].y = originVol0.y + sy * (double)dest.dimy;
        originVol[c].z = originVol0.z + sz * (double)dest.dimz;
    }

    std::vector<d4Matrix> vol2img(ic);
    std::vector<d4Vector> volOrigImg(8*ic);

    for (int im = 0; im < ic; im++)
    {
        vol2img[im] = stack.worldToImage[im] * vol2world;

        for (int c = 0; c < 8; c++)
        {
            volOrigImg[8*im + c] = vol2img[im] * originVol[c];
        }
    }

    const double s2 = sigma*sigma;

    #if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif
    FOR_ALL_VOXELS(dest)
    {
        double sum = 0.0;

        d4Vector pw(x,y,z,1.0);

        for (int im = frame0; im < ic; im++)
        {
            d4Vector pi = vol2img[im] * pw;

            for (int c = 0; c < 8; c++)
            {
                d4Vector d = pi - volOrigImg[8*im + c];

                double dx, dy;

                if (d.x < -1 || d.x > 1) dx = 0.0;
                else dx = 1.0 - std::abs(d.x);
                if (d.y < -1 || d.y > 1) dy = 0.0;
                else dy = 1.0 - std::abs(d.y);

                sum += dx*dy;
            }
        }

        dest(x,y,z) = sum;
    }
}

void BackprojectionHelper::taperEdges(Volume<RFLOAT>& vol, double rx, double ry, double rz)
{
    double mean = 0.0;

    FOR_ALL_VOXELS(vol)
    {
        mean += vol(x,y,z);
    }

    mean /= ((double)vol.dimx * (double)vol.dimy * (double)vol.dimz);

    #if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif
    FOR_ALL_VOXELS(vol)
    {
        double wx(1.0), wy(1.0), wz(1.0);

        if (x < rx) wx *= (1.0 - cos(PI * (x+1) / rx))/2.0;
        if (x >= vol.dimx - rx) wx *= (1.0 - cos(PI * (vol.dimx - x) / rx))/2.0;

        if (y < ry) wy *= (1.0 - cos(PI * (y+1) / ry))/2.0;
        if (y >= vol.dimy - ry) wy *= (1.0 - cos(PI * (vol.dimy - y) / ry))/2.0;

        if (z < rz) wz *= (1.0 - cos(PI * (z+1) / rz))/2.0;
        if (z >= vol.dimz - rz) wz *= (1.0 - cos(PI * (vol.dimz - z) / rz))/2.0;

        const double ww = wx*wy*wz;
        vol(x,y,z) = ww * vol(x,y,z) + (1.0 - ww) * mean;
    }
}
