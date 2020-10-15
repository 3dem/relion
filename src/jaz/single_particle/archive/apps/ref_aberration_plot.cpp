
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/backprojector.h>
#include <src/euler.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/jaz/single_particle/complex_io.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/resampling_helper.h>
#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/refinement_program.h>
#include <src/jaz/single_particle/parallel_ft.h>
#include <src/jaz/optics/aberration_fit.h>

#include <omp.h>

using namespace gravis;

class AberrationPlot : public RefinementProgram
{
	public:

		AberrationPlot();

			RFLOAT kmin;

			bool precomputed;
			std::string precomp;

			Image<RFLOAT> lastCos, lastSin;

		int readMoreOptions(IOParser& parser, int argc, char *argv[]);
		int _init();
		int _run();
};

AberrationPlot :: AberrationPlot()
:	RefinementProgram(true),
	precomputed(false)
{
	optReference = true;
}

int main(int argc, char *argv[])
{
	AberrationPlot tf;

	int rc0 = tf.init(argc, argv);
	if (rc0 != 0) return rc0;

	int rc1 = tf.run();
	if (rc1 != 0) return rc1;
}

int AberrationPlot::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin", "Inner freq. threshold [Angst]", "30.0"));

	precomp = parser.getOption("--precomp", "Precomputed *_sin and *_cos files from previous run (optional)", "");

	precomputed = precomp != "";

	noReference = precomputed;

	bool allGood = true;

	if (reconFn0 == "" && !precomputed)
	{
		std::cerr << "A reference map (--m) is required if no precomputed pixel-fit is available (--precomp).\n";
		allGood = false;
	}

	if (!allGood) return RELION_EXIT_FAILURE;
	else return RELION_EXIT_SUCCESS;
}

int AberrationPlot::_init()
{
	return RELION_EXIT_SUCCESS;
}

int AberrationPlot::_run()
{
	std::vector<ParFourierTransformer> fts(nr_omp_threads);

	double t0 = omp_get_wtime();

	const bool differential = false;

	if (differential)
	{
		std::vector<Image<double>>
				A(nr_omp_threads, Image<double>(sh,s)),
				b(nr_omp_threads, Image<double>(sh,s));

		const double as = (double)s * angpix;

		for (long g = minMG; g <= gc; g++)
		{
			std::stringstream stsg;
			stsg << g;

			std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

			const int pc = mdts[g].numberOfObjects();

			std::vector<Image<Complex>> pred;
			std::vector<Image<Complex>> obsF;

			pred = obsModel.predictObservations(projectors[0], mdts[g], nr_omp_threads, false, true);
			obsF = StackHelper::loadStackFS(mdts[g], imgPath, nr_omp_threads);

			#pragma omp parallel for num_threads(nr_omp_threads)
			for (long p = 0; p < pc; p++)
			{
				int t = omp_get_thread_num();

				CTF ctf0;
				ctf0.read(mdts[g], mdts[g], p);
				//ctf0.Cs = 0.0;
				ctf0.initialise();

				for (int y = 0; y < s;	y++)
				for (int x = 0; x < sh; x++)
				{
					const double xf = x;
					const double yf = y < sh? y : y - s;
					const double gamma_i = ctf0.getGamma(xf/as, yf/as);
					const double cg = cos(gamma_i);
					const double sg = sin(gamma_i);

					Complex zobs = obsF[p](y,x);
					Complex zprd = pred[p](y,x);

					double zz = zobs.real*zprd.real + zobs.imag*zprd.imag;
					double nr = zprd.norm();

					A[t](y,x) += nr*cg*cg;
					b[t](y,x) += cg*(sg*nr+zz);
				}
			}
		}

		for (int t = 1; t < nr_omp_threads; t++)
		{
			for (int y = 0; y < s;	y++)
			for (int x = 0; x < sh; x++)
			{
				A[0](y,x) += A[t](y,x);
				b[0](y,x) += b[t](y,x);
			}
		}

		Image<RFLOAT> dgamma(sh,s);

		for (int y = 0; y < s;	y++)
		for (int x = 0; x < sh; x++)
		{
			if (A[0](y,x) != 0.0)
			{
				dgamma(y,x) = b[0](y,x)/A[0](y,x);
			}
		}

		ImageLog::write(dgamma, outPath+"_dgamma");
	}
	else
	{
		Image<RFLOAT> cosPhi(sh,s), sinPhi(sh,s), phase(sh,s);

		if (precomputed)
		{
			std::cout << "using precomputed data...\n";

			cosPhi.read(precomp+"_cos.mrc");
			sinPhi.read(precomp+"_sin.mrc");
			phase.read(precomp+"_phase.mrc");

			s = cosPhi.data.ydim;
			sh = cosPhi.data.xdim;
		}
		else
		{
			std::vector<Image<double>>
				Axx(nr_omp_threads, Image<double>(sh,s)),
				Axy(nr_omp_threads, Image<double>(sh,s)),
				Ayy(nr_omp_threads, Image<double>(sh,s)),
				bx(nr_omp_threads, Image<double>(sh,s)),
				by(nr_omp_threads, Image<double>(sh,s));

			const double as = (double)s * angpix;

			for (long g = minMG; g <= gc; g++)
			{
				std::stringstream stsg;
				stsg << g;

				std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

				const int pc = mdts[g].numberOfObjects();

				std::vector<Image<Complex> > pred;
				std::vector<Image<Complex> > obsF;

				pred = obsModel.predictObservations(projectors[0], mdts[g], nr_omp_threads, false, true);
				obsF = StackHelper::loadStackFS(mdts[g], imgPath, nr_omp_threads);

				#pragma omp parallel for num_threads(nr_omp_threads)
				for (long p = 0; p < pc; p++)
				{
					int t = omp_get_thread_num();

					CTF ctf0;
					ctf0.read(mdts[g], mdts[g], p);
					//ctf0.Cs = 0.0;
					ctf0.initialise();

					for (int y = 0; y < s;	y++)
					for (int x = 0; x < sh; x++)
					{
						const double xf = x;
						const double yf = y < sh? y : y - s;
						const double gamma_i = ctf0.getGamma(xf/as, yf/as);
						const double cg = cos(gamma_i);
						const double sg = sin(gamma_i);

						Complex zobs = obsF[p](y,x);
						Complex zprd = pred[p](y,x);

						double zz = zobs.real*zprd.real + zobs.imag*zprd.imag;
						double nr = zprd.norm();

						Axx[t](y,x) += nr*sg*sg;
						Axy[t](y,x) += nr*cg*sg;
						Ayy[t](y,x) += nr*cg*cg;

						bx[t](y,x) -= zz*sg;
						by[t](y,x) -= zz*cg;
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
					Ayy[0](y,x) += Ayy[t](y,x);

					bx[0](y,x) += bx[t](y,x);
					by[0](y,x) += by[t](y,x);
				}
			}

			for (int y = 0; y < s;	y++)
			for (int x = 0; x < sh; x++)
			{
				d2Matrix A(
					Axx[0](y,x), Axy[0](y,x),
					Axy[0](y,x), Ayy[0](y,x));

				d2Vector b(bx[0](y,x), by[0](y,x));

				double det = A(0,0)*A(1,1) - A(1,0)*A(0,1);

				if (det != 0.0)
				{
					d2Matrix Ai = A;
					Ai.invert();

					d2Vector opt = Ai*b;

					cosPhi(y,x) = opt.x;
					sinPhi(y,x) = opt.y;

					phase(y,x) = std::abs(opt.x) > 0.0? atan2(opt.y, opt.x) : 0.0;
				}
			}

			ImageLog::write(cosPhi, outPath+"_cos");
			ImageLog::write(sinPhi, outPath+"_sin");
			//ImageLog::write(phase, outPath+"_phase");

			Image<RFLOAT> phaseFull(s,s);
			FftwHelper::decenterDouble2D(phase.data, phaseFull.data);
			ImageLog::write(phaseFull, outPath+"_phase");

			cosPhi.write(outPath+"_cos.mrc");
			sinPhi.write(outPath+"_sin.mrc");
			phase.write(outPath+"_phase.mrc");
		}

		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			double xx = x;
			double yy = y <= sh? y : y - s;
			double r = sqrt(xx*xx + yy*yy);

			if (r == 0 || 2.0*sh*angpix/r > kmin)
			{
				freqWeight(y,x) = 0.0;
			}
		}

		OriginalBasis fit = AberrationFit::fitBasic(phase, freqWeight, angpix);

		Image<RFLOAT> vis = AberrationFit::draw(&fit, angpix, s);
		Image<RFLOAT> visFull(s,s);
		FftwHelper::decenterDouble2D(vis.data, visFull.data);
		ImageLog::write(visFull, outPath+"_fit");


		MetaDataTable mdtAll;
		mdtAll.reserve(mdt0.numberOfObjects());

		for (long g = minMG; g <= gc; g++)
		{
			const int pc = mdts[g].numberOfObjects();

			for (long p = 0; p < pc; p++)
			{
				fit.offsetCtf(mdts[g], p);
			}

			mdtAll.append(mdts[g]);
		}

		mdtAll.write(outPath+".star");
	}

	double t1 = omp_get_wtime();

	std::cout << "elapsed: " << (t1 - t0) << "s \n";

	return RELION_EXIT_SUCCESS;
}
