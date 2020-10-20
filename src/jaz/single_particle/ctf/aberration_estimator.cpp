#include "aberration_estimator.h"
#include "tilt_helper.h"
#include "ctf_refiner.h"

#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/reference_map.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/complex_io.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/color_helper.h>

#include <src/args.h>
#include <src/ctf.h>

#include <set>

#include <omp.h>

using namespace gravis;


AberrationEstimator::AberrationEstimator()
	:	ready(false)
{}

void AberrationEstimator::read(IOParser &parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin_aberr",
		"Inner freq. threshold for symmetrical aberration estimation [Å]", "20.0"));

	std::string aberrToken = "--even_aberr_max_n";

	aberr_n_max  = textToInteger(parser.getOption(aberrToken,
		"Maximum degree of Zernike polynomials used to fit even (i.e. symmetrical) aberrations", "4"));

	xring0 = textToDouble(parser.getOption("--xr0_a", "Exclusion ring start [Å]", "-1"));
	xring1 = textToDouble(parser.getOption("--xr1_a", "Exclusion ring end [Å]", "-1"));
}

void AberrationEstimator::init(
		int verb, int nr_omp_threads,
		bool debug, bool diag, std::string outPath,
		ReferenceMap* reference, ObservationModel* obsModel)
{
	this->verb = verb;
	this->nr_omp_threads = nr_omp_threads;

	this->debug = debug;
	this->diag = diag;
	this->outPath = outPath;

	this->reference = reference;
	this->obsModel = obsModel;

	angpix = obsModel->getPixelSizes();
	obsModel->getBoxSizes(s, sh);

	ready = true;
}

void AberrationEstimator::processMicrograph(
		long g, MetaDataTable& mdt,
		const std::vector<Image<Complex>>& obs,
		const std::vector<Image<Complex>>& pred)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: AberrationEstimator::processMicrograph: AberrationEstimator not initialized.");
	}

	std::vector<std::pair<int, std::vector<int>>> particlesByOpticsGroup
			= obsModel->splitParticlesByOpticsGroup(mdt);

	for (int pog = 0; pog < particlesByOpticsGroup.size(); pog++)
	{
		const int og = particlesByOpticsGroup[pog].first;
		const std::vector<int>& partIndices = particlesByOpticsGroup[pog].second;

		// TODO: SHWS 29mar2018: when data is CTF-premultiplied: do we need to change below??
		if (obsModel->getCtfPremultiplied(og))
			std::cerr << "TODO: check aberration estimation with CTF-premultiplied data!!" << std::endl;

		const int pc = partIndices.size();

		std::vector<Image<RFLOAT>>
			Axx(nr_omp_threads, Image<RFLOAT>(sh[og],s[og])),
			Axy(nr_omp_threads, Image<RFLOAT>(sh[og],s[og])),
			Ayy(nr_omp_threads, Image<RFLOAT>(sh[og],s[og])),
			bx(nr_omp_threads, Image<RFLOAT>(sh[og],s[og])),
			by(nr_omp_threads, Image<RFLOAT>(sh[og],s[og]));

		const double as = (double)s[og] * angpix[og];

		#pragma omp parallel for num_threads(nr_omp_threads)
		for (long pp = 0; pp < pc; pp++)
		{
			const int p = partIndices[pp];

			CTF ctf;
			ctf.readByGroup(mdt, obsModel, p);

			int t = omp_get_thread_num();

			for (int y = 0; y < s[og];  y++)
			for (int x = 0; x < sh[og]; x++)
			{
				const double xf = x;
				const double yf = y < sh[og]? y : y - s[og];

				const double gamma_i = ctf.getLowOrderGamma(xf/as, yf/as);
				const double cg = cos(gamma_i);
				const double sg = sin(gamma_i);

				Complex zobs = obs[p](y,x);
				Complex zprd = pred[p](y,x);

				double zz = zobs.real * zprd.real + zobs.imag * zprd.imag;
				double nr = zprd.norm();

				Axx[t](y,x) += nr * sg * sg;
				Axy[t](y,x) += nr * cg * sg;
				Ayy[t](y,x) += nr * cg * cg;

				bx[t](y,x) -= zz * sg;
				by[t](y,x) -= zz * cg;
			}
		}

		// Combine the accumulated weights from all threads for this subset

		Image<RFLOAT>
			AxxSum(sh[og],s[og]), AxySum(sh[og],s[og]), AyySum(sh[og],s[og]),
			bxSum(sh[og],s[og]), bySum(sh[og],s[og]);

		for (int threadnum = 0; threadnum < nr_omp_threads; threadnum++)
		{
			ImageOp::linearCombination(AxxSum, Axx[threadnum], 1.0, 1.0, AxxSum);
			ImageOp::linearCombination(AxySum, Axy[threadnum], 1.0, 1.0, AxySum);
			ImageOp::linearCombination(AyySum, Ayy[threadnum], 1.0, 1.0, AyySum);

			ImageOp::linearCombination(bxSum, bx[threadnum], 1.0, 1.0, bxSum);
			ImageOp::linearCombination(bySum, by[threadnum], 1.0, 1.0, bySum);
		}

		// Write out the intermediate results per-micrograph:

		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);

		std::stringstream sts;
		sts << (og+1);

		AxxSum.write(outRoot+"_aberr-Axx_optics-group_" + sts.str() + ".mrc");
		AxySum.write(outRoot+"_aberr-Axy_optics-group_" + sts.str() + ".mrc");
		AyySum.write(outRoot+"_aberr-Ayy_optics-group_" + sts.str() + ".mrc");

		bxSum.write(outRoot+"_aberr-bx_optics-group_" + sts.str() + ".mrc");
		bySum.write(outRoot+"_aberr-by_optics-group_" + sts.str() + ".mrc");
	}
}

void AberrationEstimator::parametricFit(
		const std::vector<MetaDataTable>& mdts,
		MetaDataTable& optOut, std::vector <FileName> &fn_eps)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: AberrationEstimator::parametricFit: AberrationEstimator not initialized.");
	}

	if (verb > 0)
	{
		std::cout << " + Fitting symmetrical aberrations ..." << std::endl;
	}

	const int gc = mdts.size();
	const int ogc = obsModel->numberOfOpticsGroups();

	std::vector<bool> groupUsed(ogc,false);

	#pragma omp parallel for num_threads(nr_omp_threads)
	for (int og = 0; og < ogc; og++)
	{
		std::stringstream sts;
		sts << og+1;
		std::string ogstr = sts.str();

		Image<RFLOAT>
			AxxSum(sh[og],s[og]), AxySum(sh[og],s[og]), AyySum(sh[og],s[og]),
			bxSum(sh[og],s[og]), bySum(sh[og],s[og]);

		for (long g = 0; g < gc; g++)
		{
			std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);

			if (   exists(outRoot+"_aberr-Axx_optics-group_" + ogstr + ".mrc")
				&& exists(outRoot+"_aberr-Axy_optics-group_" + ogstr + ".mrc")
				&& exists(outRoot+"_aberr-Ayy_optics-group_" + ogstr + ".mrc")
				&& exists(outRoot+"_aberr-bx_optics-group_" + ogstr + ".mrc")
				&& exists(outRoot+"_aberr-by_optics-group_" + ogstr + ".mrc"))
			{

				Image<RFLOAT>
					Axx(sh[og],s[og]), Axy(sh[og],s[og]), Ayy(sh[og],s[og]),
					bx(sh[og],s[og]), by(sh[og],s[og]);

				Axx.read(outRoot+"_aberr-Axx_optics-group_" + ogstr + ".mrc");
				Axy.read(outRoot+"_aberr-Axy_optics-group_" + ogstr + ".mrc");
				Ayy.read(outRoot+"_aberr-Ayy_optics-group_" + ogstr + ".mrc");

				bx.read(outRoot+"_aberr-bx_optics-group_" + ogstr + ".mrc");
				by.read(outRoot+"_aberr-by_optics-group_" + ogstr + ".mrc");

				AxxSum() += Axx();
				AxySum() += Axy();
				AyySum() += Ayy();

				bxSum() += bx();
				bySum() += by();

				groupUsed[og] = true;
			}
		}

		if (!groupUsed[og])
		{
			continue;
		}

		Image<RFLOAT> wgh0(sh[og],s[og]), wgh(sh[og],s[og]), phase(sh[og],s[og]);
		Image<Complex> optXY(sh[og],s[og]);

		wgh0 = reference->getHollowWeight(kmin, s[og], angpix[og]);

		for (int y = 0; y < s[og];  y++)
		for (int x = 0; x < sh[og]; x++)
		{
			d2Matrix A(
				AxxSum(y,x), AxySum(y,x),
				AxySum(y,x), AyySum(y,x));

			d2Vector b(bxSum(y,x), bySum(y,x));

			double det = A(0,0) * A(1,1) - A(1,0) * A(0,1);

			if (det != 0.0)
			{
				d2Matrix Ai = A;
				Ai.invert();

				d2Vector opt = Ai * b;

				optXY(y,x) = Complex(opt.x, opt.y);
				phase(y,x) = std::abs(opt.x) > 0.0? atan2(opt.y, opt.x) : 0.0;
				wgh(y,x) = wgh0(y,x) * sqrt(std::abs(det));
			}
			else
			{
				optXY(y,x) = 0.0;
				phase(y,x) = 0.0;
				wgh0(y,x) = 0.0;
				wgh(y,x) = 0.0;
			}
		}

		if (xring1 > 0.0)
		{
			for (int y = 0; y < s[og]; y++)
			for (int x = 0; x < sh[og]; x++)
			{
				double xx = x;
				double yy = y < sh[og]? y : y - s[og];
				double rp = sqrt(xx*xx + yy*yy);
				double ra = s[og] * angpix[og] / rp;

				if (ra > xring0 && ra <= xring1)
				{
					wgh0(y,x) = 0.0;
					wgh(y,x) = 0.0;
				}
			}
		}

		if (debug)
		{
			Image<RFLOAT> full;
			FftwHelper::decenterDouble2D(wgh(), full());
			ImageLog::write(full, outPath + "aberr_weight-full_optics-group_"+ogstr);
		}

		std::vector<Image<RFLOAT> > imgs_for_eps;
		std::vector<double> scales;
		std::vector<std::string> labels;

		Image<RFLOAT> fit, phaseFull, fitFull;
		FftwHelper::decenterDouble2D(phase.data, phaseFull.data);
		ImageLog::write(phaseFull, outPath + "aberr_delta-phase_per-pixel_optics-group_"+ogstr);

		imgs_for_eps.push_back(phaseFull);
		scales.push_back(1.);
		labels.push_back("Symm. obs [-1, 1] "+obsModel->getGroupName(og));
		imgs_for_eps.push_back(phaseFull);
		scales.push_back(PI);
		labels.push_back("Symm. obs [-pi, pi] "+obsModel->getGroupName(og));


		{
			std::vector<double> Zernike_coeffs = TiltHelper::fitEvenZernike(
				phase, wgh, angpix[og], obsModel->getMagMatrix(og), aberr_n_max, &fit);

			FftwHelper::decenterDouble2D(fit.data, fitFull.data);

			std::stringstream sts;
			sts << aberr_n_max;

			ImageLog::write(fitFull, outPath + "aberr_delta-phase_lin-fit_optics-group_"
							+ogstr+"_N-"+sts.str());

			{
				Image<RFLOAT> residual;
				residual.data = phaseFull.data - fitFull.data;

				ImageLog::write(residual, outPath + "aberr_delta-phase_lin-fit_optics-group_"
								+ogstr+"_N-"+sts.str()+"_residual");
			}

			std::vector<double> Zernike_coeffs_opt = TiltHelper::optimiseEvenZernike(
				optXY, wgh0, AxxSum, AxySum, AyySum,
				angpix[og], obsModel->getMagMatrix(og),
				aberr_n_max, Zernike_coeffs, &fit);

			FftwHelper::decenterDouble2D(fit.data, fitFull.data);

			ImageLog::write(fitFull, outPath + "aberr_delta-phase_iter-fit_optics-group_"
							+ogstr+"_N-"+sts.str());

			imgs_for_eps.push_back(fitFull);
			scales.push_back(1.);
			labels.push_back("Symm. (N="+sts.str()+") fit [-1, 1] "+obsModel->getGroupName(og));
			imgs_for_eps.push_back(fitFull);
			scales.push_back(PI);
			labels.push_back("Symm. (N="+sts.str()+") fit [-pi, pi] "+obsModel->getGroupName(og));


			// extract Q0, Cs, defocus and astigmatism?
			#pragma omp critical
			{
				optOut.setValue(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, Zernike_coeffs_opt, og);
			}
		}

		FileName fn_root = outPath + "symmetric_aberrations_optics-group_"+ ogstr;
		ColorHelper::writeSignedToEPS(fn_root, 2, imgs_for_eps, scales, labels);
		fn_eps.push_back(fn_root+".eps");

	}

}

bool AberrationEstimator::isFinished(const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: AberrationEstimator::isFinished: AberrationEstimator not initialized.");
	}

	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);

	bool allDone = true;

	std::vector<int> ogs = obsModel->getOptGroupsPresent_oneBased(mdt);

	for (int i = 0; i < ogs.size(); i++)
	{
		const int og = ogs[i];

		std::stringstream sts;
		sts << og;

		if (   !exists(outRoot+"_aberr-Axx_optics-group_" + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-Axy_optics-group_" + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-Ayy_optics-group_" + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-bx_optics-group_"  + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-by_optics-group_"  + sts.str() + ".mrc"))
		{
			allDone = false;
			break;
		}
	}

	return allDone;
}
