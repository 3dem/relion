#ifndef ABERRATION_PLOT_H
#define ABERRATION_PLOT_H

#include <string>

#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/backprojector.h>

#include <src/jaz/obs_model.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/gravis/t3Matrix.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;

	parser.setCommandLine(argc, argv);
	parser.addSection("General options");

	std::string starFn = parser.getOption("--i", "Input STAR file");
	std::string outPath = parser.getOption("--out", "Output path");
	RFLOAT angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "0.0"));

	std::string imgPath = parser.getOption("--img", "Path to images", "");

	int nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of OMP threads", "1"));
	long minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
	long maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index", "-1"));

	if (parser.checkForErrors()) return 1;

	std::cout << "reading " << starFn << "...\n";

	MetaDataTable mdt0;
	mdt0.read(starFn);

	if (angpix <= 0.0)
	{
		RFLOAT mag, dstep;
		mdt0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
		mdt0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
		angpix = 10000 * dstep / mag;

		std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << "\n";
	}

	std::vector<MetaDataTable> mdts;
	mdts = StackHelper::splitByStack(&mdt0);


	std::string name, imgName;
	mdts[0].getValue(EMDL_IMAGE_NAME, imgName, 0);
	name = imgName.substr(imgName.find("@")+1);

	std::string finName;

	if (imgPath == "")
	{
		finName = name;
	}
	else
	{
		finName = imgPath + "/" + imgName.substr(imgName.find_last_of("/")+1);
	}

	Image<RFLOAT> stack0;
	stack0.read(finName, false);

	const int s = stack0.data.xdim;
	const int sh = s/2 + 1;


	const long gc = maxMG >= 0? maxMG : mdts.size()-1;
	const long g0 = minMG;

	std::cout << "mg range: " << g0 << ".." << gc << "\n";

	std::vector<ParFourierTransformer> fts(nr_omp_threads);

	std::vector<Image<double>>
			Axx(nr_omp_threads, Image<double>(sh,s)),
			Axy(nr_omp_threads, Image<double>(sh,s)),
			Axz(nr_omp_threads, Image<double>(sh,s)),
			Ayy(nr_omp_threads, Image<double>(sh,s)),
			Ayz(nr_omp_threads, Image<double>(sh,s)),
			Azz(nr_omp_threads, Image<double>(sh,s)),
			bx(nr_omp_threads, Image<double>(sh,s)),
			by(nr_omp_threads, Image<double>(sh,s)),
			bz(nr_omp_threads, Image<double>(sh,s));


	const double as = (double)s * angpix;

	double t0 = omp_get_wtime();

	for (long g = minMG; g <= gc; g++)
	{
		std::stringstream stsg;
		stsg << g;

		std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

		const int pc = mdts[g].numberOfObjects();

		std::vector<Image<Complex> > obsF;

		obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

		#pragma omp parallel for num_threads(nr_omp_threads)
		for (long p = 0; p < pc; p++)
		{
			int t = omp_get_thread_num();

			CTF ctf0;
			ctf0.read(mdts[g], mdts[g], p);
			ctf0.Cs = 0.0;
			ctf0.initialise();

			for (int y = 0; y < s;	y++)
			for (int x = 0; x < sh; x++)
			{
				const double xf = x;
				const double yf = y < sh? y : y - s;
				const double gamma = ctf0.getGamma(xf/as, yf/as);
				const double cg = cos(gamma);
				const double sg = sin(gamma);

				double cx = cg*cg;
				double cy = 2.0*sg*cg;
				double cz = sg*sg;

				double zz = obsF[p](y,x).norm();

				Axx[t](y,x) += cx*cx;
				Axy[t](y,x) += cx*cy;
				Axz[t](y,x) += cx*cz;
				Ayy[t](y,x) += cy*cy;
				Ayz[t](y,x) += cy*cz;
				Azz[t](y,x) += cz*cz;

				bx[t](y,x) += zz*cx;
				by[t](y,x) += zz*cy;
				bz[t](y,x) += zz*cz;
			}
		}
	}

	for (int t = 1; t < nr_omp_threads; t++)
	{
		for (int y = 0; y < s;	y++)
		for (int x = 0; x < sh; x++)
		{
			Axx[0](y,x) += Axx[t](y,x);
			Axy[0](y,x) += Axy[t](y,x);
			Axz[0](y,x) += Axz[t](y,x);
			Ayy[0](y,x) += Ayy[t](y,x);
			Ayz[0](y,x) += Ayz[t](y,x);
			Azz[0](y,x) += Azz[t](y,x);

			bx[0](y,x) += bx[t](y,x);
			by[0](y,x) += by[t](y,x);
			bz[0](y,x) += bz[t](y,x);
		}
	}

	Image<RFLOAT> xx(sh,s), xy(sh,s), xz(sh,s);

	for (int y = 0; y < s;	y++)
	for (int x = 0; x < sh; x++)
	{
		d3Matrix A(
		           Axx[0](y,x), Axy[0](y,x), Axz[0](y,x),
		           Axy[0](y,x), Ayy[0](y,x), Ayz[0](y,x),
		           Axz[0](y,x), Ayz[0](y,x), Azz[0](y,x));

		d3Vector b(bx[0](y,x), by[0](y,x), bz[0](y,x));

		d3Matrix Ai = A;
		Ai.invert();

		d3Vector opt = Ai*b;

		xx(y,x) = opt.x;
		xy(y,x) = opt.y;
		xz(y,x) = opt.z;
	}

	ImageLog::write(xx, outPath+"_sin2");
	ImageLog::write(xy, outPath+"_sincos");
	ImageLog::write(xz, outPath+"_cos2");

	return 0;
}

#endif
