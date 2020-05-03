#include "color_helper.h"
#include <src/macros.h>
#include <cmath>

#ifdef HAVE_PNG
#include <src/jaz/gravis/tImage.h>
#endif

using namespace gravis;

dRGB ColorHelper::signedToRedBlue(double d, double scale, double rbFract)
{		
	const double d_rb = d / (scale * rbFract);
	const double d_g = (std::abs(d)/scale - rbFract) / (1.0 - rbFract);
	
	return dRGB(
		std::min(1.0, std::max(0.0,  d_rb)), 
		std::min(1.0, std::max(0.0,  d_g)), 
		std::min(1.0, std::max(0.0, -d_rb)) );
}

void ColorHelper::writeAngleToPNG(const Image<RFLOAT> &img, std::string filename)
{
	writeSignedToPNG(img, filename+"_[-pi,+pi]", PI);
	writeSignedToPNG(img, filename+"_[-1,+1]", 1.0);
}

void ColorHelper::writeSignedToPNG(const Image<RFLOAT> &img, std::string filename, double scale)
{
	#ifdef HAVE_PNG
	{
		tImage<dRGB> pngOut(img.data.xdim, img.data.ydim);
		pngOut.fill(dRGB(0.f));
		
		for (int y = 0; y < img.data.ydim; y++)
		for (int x = 0; x < img.data.xdim; x++)
		{
			double c = img(y,x);
			pngOut(x,y) = signedToRedBlue(c, scale);
		}
		
		pngOut.writePNG(filename+".png");
	}
	#endif
}

void ColorHelper::writeSignedToEPS(std::string filename, int col, const std::vector<Image<RFLOAT> > &imgs,
		const std::vector<double> &scales, const std::vector<std::string> &labels)
{
	// Check all images have the same size
	int xdim = imgs[0].data.xdim;
	int ydim = imgs[0].data.ydim;
	int nimgs= imgs.size();
	for (int i = 1; i < imgs.size(); i++)
	{
		if (imgs[i].data.xdim != xdim || imgs[i].data.ydim != ydim)
			REPORT_ERROR(" ERROR: combining images with different sizes into one EPS...");
	}

	std::ofstream outputFile;
	FileName fn_out = filename + ".eps";
	outputFile.open(fn_out.c_str());

	int delta = 15;
	int row = CEIL(nimgs/(RFLOAT)col);
	int width = col * xdim + (col-1)*delta;
	int height = row * (ydim + delta);

	// Rescale to maximum one A4: 595 x 842, or one letter: 612 x 792
	RFLOAT width_ratio = width / 595.;
	RFLOAT height_ratio = height / 792.;
	RFLOAT max_ratio = XMIPP_MAX(width_ratio, height_ratio);
	RFLOAT rescale = 1.;
	if (max_ratio > 1.)
	{
		rescale = max_ratio;
	}

	// header
	outputFile << "%!PS-Adobe-2.0 EPSF-1.2" << "\n";
	outputFile << "%%BoundingBox: 0 0 " << ROUND(width/rescale)<< " " <<  ROUND(height/rescale) << "\n";
	outputFile << "%%Pages: 1" << "\n";
	outputFile << "%%EndComments" << "\n";
	outputFile << "/Times-Roman findfont\n";
	outputFile << ROUND(10/rescale) << " scalefont\n";
	outputFile << "setfont\n";

	// First put all the labels (without scale argument!)
	int xcoord, ycoord, xpos, ypos;
	for (int i = 0; i < imgs.size(); i++)
	{
		xpos = i%col;
		ypos = (row - 1) - i/col;
		xcoord = xpos * ROUND((xdim+delta)/rescale);
		ycoord = ypos * ROUND((ydim+delta)/rescale) + ROUND(ydim/rescale);

		// Print the label
		outputFile << "newpath\n";
		outputFile << (int)(xcoord) << " " << (int)(ycoord + ROUND(5/rescale)) << " moveto\n";
		outputFile << "(" << labels[i] << ") show\n";
	}

	// one scale statement only!
	outputFile << ROUND(xdim/rescale) << " " << ROUND(ydim/rescale) << "  scale\n";
	for (int i = 0; i < imgs.size(); i++)
	{
		xpos = i%col;
		ypos = (row - 1) - i/col;
		xcoord = xpos * (xdim + delta);
		ycoord = ypos * (ydim+delta) + ydim;

		// The actual image
		// Note that the number of elements in a string or array literal should be less than 64 K.
		// Otherwise, earlier versions of Ghostscript and Preview in MacOS fails.
		// Ref: https://stackoverflow.com/questions/7595532/postcript-maximum-array-size
		//      http://paulbourke.net/dataformats/postscript/
		outputFile << xdim << " " << ydim <<" 8 [" << xdim << " 0 0 -" << ydim << " -"<<xcoord << " " << ycoord<<"]\n";
		outputFile << "{currentfile " << (xdim * 3) << " string readhexstring pop} bind\n";
		outputFile << "false 3 colorimage\n";

		long ii=0;
		for (int y = 0; y < ydim; y++)
		for (int x = 0; x < xdim; x++)
		{
			double c = imgs[i](y,x);
			dRGB myRGB = signedToRedBlue(c, scales[i]);
			outputFile << std::hex << std::setfill('0') << std::setw(2) <<(int)(myRGB.r*255);
			outputFile << std::hex << std::setfill('0') << std::setw(2) <<(int)(myRGB.g*255);
			outputFile << std::hex << std::setfill('0') << std::setw(2) <<(int)(myRGB.b*255);
			ii++;
			if (ii > 6)
			{
				ii=0;
				outputFile << "\n";
			}
		}
		if (ii!=0) outputFile << "\n";
		outputFile << std::dec;
		outputFile << "\n";
	}

	outputFile << "%%EOF\n";

	outputFile.close();
}

