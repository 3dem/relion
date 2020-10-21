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

#ifndef BACKPROJECTON_HELPER_H
#define BACKPROJECTON_HELPER_H

#include <src/image.h>
#include <src/jaz/single_particle/volume.h>
#include "tomo_stack.h"
#include <src/jaz/gravis/t3Vector.h>
#include <string>

class BackprojectionHelper
{
    public:

        enum InterpolationType {Linear, Cubic};


        static void backprojectRaw(
						const Image<RFLOAT>& stack, std::string tiltAngles,
						Volume<RFLOAT>& dest, Volume<unsigned char>& maskDest,
						gravis::d3Vector origin, double spacing = 1.0, int frames = -1);

        static void backprojectRaw(
						const TomoStack& stack,
						Volume<RFLOAT>& dest, Volume<RFLOAT>& maskDest,
						gravis::d3Vector origin, double spacing = 1.0, 
						InterpolationType interpolation = Linear,
						double taperX = 20, double taperY = 20, 
						double wMin = 3.0, int frame0 = 0, int frames = -1);

        static void backprojectExactWeights(const TomoStack& stack,
                                Volume<RFLOAT>& dest,
                                gravis::d3Vector origin, double spacing = 1.0,
                                double taperX = 20, double taperY = 20, double taperZ = 20, double wMin = 3.0,
                                int frame0 = 0, int frames = -1);

        static void backprojectExactWeightsFreq(const TomoStack& stack,
                                Image<Complex>& dest,
                                Volume<RFLOAT>& weight,
                                gravis::d3Vector origin, double spacing = 1.0,
                                double taperX = 20, double taperY = 20, double taperZ = 20, double wMin = 3.0,
                                int frame0 = 0, int frames = -1);

        static void backprojectDots(const TomoStack& stack,
                                Volume<RFLOAT>& dest, gravis::d3Vector origin, double spacing = 1.0,
                                double taperX = 20, double taperY = 20, double taperZ = 20,
                                int frame0 = 0, int frames = -1);

        static void backprojectDotsFS(const TomoStack& stack,
                                Image<Complex>& dest, gravis::d3Vector origin, double spacing = 1.0,
                                double taperX = 20, double taperY = 20, double taperZ = 20,
                                int frame0 = 0, int frames = -1);

        static void backprojectDotsSeparately(const TomoStack& stack,
                                Volume<RFLOAT>& dest, gravis::d3Vector origin, double spacing = 1.0,
                                double taperX = 20, double taperY = 20, double taperZ = 20,
                                int frame0 = 0, int frames = -1);

        static void backprojectOriginDot(const TomoStack& stack,
                                Volume<RFLOAT>& dest, double sigma,
                                gravis::d3Vector origin, double spacing = 1.0,
                                int frame0 = 0, int frames = -1);

        static void taperEdges(Volume<RFLOAT>& vol, double rx, double ry, double rz);

};

#endif
