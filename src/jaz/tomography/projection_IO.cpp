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

#include "projection_IO.h"
#include <fstream>
#include <cmath>
#include <src/error.h>
#include <src/jaz/gravis/t2Vector.h>
#include <string.h>

using namespace gravis;

#define PROJ_VERSION 2

std::vector<gravis::d4Matrix> ProjectionIO::loadTiltProjections(
        std::string tiltFile, 
		double centerX, double centerY, 
		bool flipYZ, bool flipY)
{
    std::ifstream anglesFile(tiltFile.c_str());

    if (!anglesFile.is_open())
    {
        REPORT_ERROR("failed to open " + tiltFile + '\n');
    }

    std::vector<double> angles;
    std::vector<d4Matrix> vol2img;

    const double deg2rad = 3.14159265358979323846/180.0;
	
	const int inX = 0;
	const int inY = flipYZ? 2 : 1;
	const int inZ = flipYZ? 1 : 2;
	
	const int outX = 0;
	const int outY = flipYZ? 1 : 2;
	const int outZ = flipYZ? 2 : 1;
	
	d4Matrix w2i0;
    w2i0(inX,3) = -centerX;
    w2i0(inY,3) = -centerY;
	
	const double zSign = flipY? -1.0 : 1.0;
	
    while (anglesFile.good())
    {
        double a;
        anglesFile >> a;
        a *= deg2rad;

        angles.push_back(a);

        d4Matrix w2i;
		
		w2i(outX, inX) =  cos(a);
		w2i(outX, inY) =  0;
        w2i(outX, inZ) =  sin(a) * zSign;
		
		w2i(outY, inX) =  0;
		w2i(outY, inY) =  1;
        w2i(outY, inZ) =  0;
		
		w2i(outZ, inX) = -sin(a);
		w2i(outZ, inY) =  0;
        w2i(outZ, inZ) =  cos(a) * zSign;
		
        w2i(0,3) =  centerX;
		w2i(1,3) =  centerY;
		
		d4Matrix A = w2i * w2i0;

        vol2img.push_back(A);
    }

    return vol2img;
}

std::vector<d4Matrix> ProjectionIO::loadAffineTransforms(
		std::string xformFile, double cx, double cy, bool square_result)
{
    std::ifstream file(xformFile.c_str());

    if (!file.is_open())
    {
        REPORT_ERROR("failed to open " + xformFile + '\n');
    }

    std::vector<d4Matrix> xforms;

    d4Matrix P, Q;
    P.loadIdentity();
    Q.loadIdentity();

    P(0,3) = -cx;
    P(1,3) = -cy;
    Q(0,3) =  cx;

    if (square_result)
    {
        Q(1,3) =  cx;
    }
    else
    {
        Q(1,3) =  cy;
    }

    char text[4096];

    while (file.good())
    {
        file.getline(text, 4096);

        if (strlen(text) < 11) break;

        std::stringstream line(text);

        d4Matrix A;
        A.loadIdentity();

        line >> A(0,0);
        line >> A(0,1);
        line >> A(1,0);
        line >> A(1,1);
        line >> A(0,3);
        line >> A(1,3);

        xforms.push_back(Q*A*P);
    }

	return xforms;
}

void ProjectionIO::write(std::vector<d4Matrix> projs, int w, int h, int d, std::string filename)
{
	const int fc = projs.size();
	
	std::ofstream outFile(filename);

    outFile << "v " << PROJ_VERSION << "\n\n";
    outFile << fc << "\n\n";
	
	
	for (int f = 0; f < fc; f++)
    {
		outFile << f << "\n";	
				   
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				outFile << std::setprecision(12);
				outFile << std::setw(12);
				outFile << projs[f](j,k);
				
				if (k < 3) 
				{
					outFile << " ";		
				}
			}
			
			outFile << "\n";
		}
		
		outFile << "\n";
	}

    outFile << w << " " << h << " " << d << "\n";
}

std::vector<d4Matrix> ProjectionIO::read(std::string filename, int& w, int& h, int& d)
{
	std::ifstream inFile(filename);
	
	if (!inFile.good())
	{
		REPORT_ERROR("Unable to read " + filename);
	}

    int version;
    std::string v;
    inFile >> v;
    inFile >> version;

    if (v != "v" || version != PROJ_VERSION)
    {
        std::cerr << "\nERROR: The file '" << filename << "' does not have the right version number (" << PROJ_VERSION << ").\n" << std::endl;
		std::exit(RELION_EXIT_FAILURE);
    }
	
	int fc = -1;	
	inFile >> fc;
	
	std::vector<d4Matrix> out(fc);
		
	for (int ff = 0; ff < fc; ff++)
    {
		int f = -1;
		inFile >> f;
						   
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				inFile >> out[f](j,k);
			}
		}
	}
	
	inFile >> w;
	inFile >> h;
	inFile >> d;

	return out;
}

void ProjectionIO::read(
		std::string projFn, int fc, std::string stackFn, 
		std::vector<d4Matrix>& matrices_out, 
		d3Vector& centre_out,
		int& w0_out,
		int& h0_out,
		int& d0_out)
{
	matrices_out = ProjectionIO::read(projFn, w0_out, h0_out, d0_out);
	
	if (fc > 0 && matrices_out.size() != fc)
	{
		REPORT_ERROR_STR(stackFn << " contains " << fc << " frames, while " << projFn
						 << " contains " << matrices_out.size());
	}
	
	centre_out = d3Vector(w0_out/2.0, h0_out/2.0, d0_out/2.0);
}

d3Vector ProjectionIO::findSampleNormal(const std::vector<d4Matrix> &proj)
{
	d3Vector sum(0.0, 0.0, 0.0);
	
	for (int f = 0; f < proj.size(); f++)
	{
		sum += d3Vector(proj[f](2,0), proj[f](2,1), proj[f](2,2));
	}
	
	d3Vector out(0.0, 0.0, 0.0);
	
	if (std::abs(sum[0]) > std::abs(sum[1]))
	{
		if (std::abs(sum[0]) > std::abs(sum[2]))
		{
			out[0] = 1.0;
		}
		else
		{
			out[2] = 1.0;
		}
	}
	else
	{
		if (std::abs(sum[1]) > std::abs(sum[2]))
		{
			out[1] = 1.0;
		}
		else
		{
			out[2] = 1.0;
		}
	}
	
	return out;
}

double ProjectionIO::getDefocusOffset(
	double x, double y, 
	int w, int h, 
	const d4Matrix& P, 
	double pixSize, 
	const d3Vector& sampleNormal)
{
	d2Vector dxy = d2Vector(x - w/2.0, y - h/2.0);
	d4Vector r4 = P * d4Vector(sampleNormal.x, sampleNormal.y, sampleNormal.z, 0.0);
	
	d2Vector r(r4.x / r4.z, r4.y / r4.z);
	
	return -pixSize * r.dot(dxy);
}
