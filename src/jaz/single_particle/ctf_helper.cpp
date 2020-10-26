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

#include <src/jaz/single_particle/ctf_helper.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>


std::vector<CTF> CtfHelper :: loadCtffind4(
		std::string path, int imageCount, 
		double voltage, double Cs, double Q0, double Bfac, double scale)
{
    /*
     example:
    # Output from CTFFind version 4.1.5, run on 2017-03-30 15:12:45
    # Input file: /beegfs/zivanov/tomograms/ts_05/frames/05_f32.mrc ; Number of micrographs: 1
    # Pixel size: 1.000 Angstroms ; acceleration voltage: 300.0 keV ; spherical aberration: 2.70 mm ; amplitude contrast: 0.07
    # Box size: 512 pixels ; min. res.: 30.0 Angstroms ; max. res.: 5.0 Angstroms ; min. def.: 5000.0 um; max. def. 50000.0 um
    # Columns: #1 - micrograph number; #2 - defocus 1 [Angstroms]; #3 - defocus 2; #4 - azimuth of astigmatism; #5 - additional phase shift [radians]; #6 - cross correlation; #7 - spacing (in Angstroms) up to which CTF rings were fit successfully
    1.000000 10295.926758 10012.275391 -38.856349 0.000000 0.030650 5.279412
    */
	
	std::vector<CTF> ctfs(imageCount);
	
    size_t ast = path.find_first_of('*');
	
    if (ast == std::string::npos)
    {
		std::ifstream file(path);
		int currImg = 0;
		
		char text[4096];

		while (file.getline(text, 4096))
		{
			if (text[0] == '#') continue;
			
			std::stringstream line(text);

			ctfs[currImg] = setFromFile(line, voltage, Cs, Q0, Bfac, scale);
			currImg++;
			
			if (currImg >= imageCount)
			{
				break;
			}
		}
		
		if (currImg < imageCount)
		{
			REPORT_ERROR_STR("Insufficient number of CTFs found in " << path << ".\n"
							 << imageCount << " requested, " << currImg << " found.\n");
		}
    }
	else
	{
		std::string fnBase, fnEnd;
		fnBase = path.substr(0, ast);
		fnEnd = path.substr(ast+1);
		
		for (int i = 0; i < imageCount; i++)
		{
			std::stringstream sts;
			sts << i;
			std::string fnm;
			sts >> fnm;
	
			std::string fn = fnBase + fnm + fnEnd;
			std::ifstream file(fn.c_str());
	
			if (!file.is_open())
			{
				REPORT_ERROR("failed to open " + fn + '\n');
			}
			
			char text[4096];
	
			while (file.getline(text, 4096))
			{
				if (text[0] == '#') continue;
	
				std::stringstream line(text);
	
				ctfs[i] = setFromFile(line, voltage, Cs, Q0, Bfac, scale);
			}
		}
	}
	
	return ctfs;
}

CTF CtfHelper::setFromFile(std::stringstream& line,
						   double voltage, double Cs, double Q0, double Bfac, double scale)
{
	/*
	#1 - micrograph number;
	#2 - defocus 1 [Angstroms];
	#3 - defocus 2;
	#4 - azimuth of astigmatism;
	#5 - additional phase shift [radians];
	#6 - cross correlation;
	#7 - spacing (in Angstroms) up to which CTF rings were fit successfully
	*/

	double imgNumber, defocus1, defocus2, azimuth, phaseShift, crossCorr, bestBefore;

	line >> imgNumber;
	line >> defocus1;
	line >> defocus2;
	line >> azimuth;
	line >> phaseShift;
	line >> crossCorr;
	line >> bestBefore;

	CTF ctf;
	ctf.setValues(defocus1, defocus2, azimuth, voltage, Cs, Q0, Bfac, scale, phaseShift);
	
	return ctf;
}
