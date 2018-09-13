#include "bfactor_refiner.h"
#include "ctf_refiner.h"

#include <omp.h>

#include <src/ctf.h>
#include <src/time.h>

#include <src/jaz/obs_model.h>
#include <src/jaz/reference_map.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/img_proc/image_op.h>
#include <src/jaz/img_proc/radial_avg.h>


using namespace gravis;

BFactorRefiner::BFactorRefiner()
:	ready(false)
{}

void BFactorRefiner::read(IOParser &parser, int argc, char *argv[])
{
	/*perMicrograph = parser.checkOption("--bfac_per_mg", 
		"Estimate B-factors per micrograph, instead of per particle");*/
	
	min_B = textToDouble(parser.getOption("--bfac_min_B",
		"Minimal allowed B-factor", "-20"));
	
	max_B = textToDouble(parser.getOption("--bfac_max_B",
		"Maximal allowed B-factor", "200"));
	
	min_scale = textToDouble(parser.getOption("--bfac_min_scale",
		"Minimal allowed scale-factor (essential for outlier rejection)", "0.2"));
	
	kmin = textToDouble(parser.getOption("--kmin_bfac",
		"Inner freq. threshold for B-factor estimation [Angst]", "30.0"));
}

void BFactorRefiner::init(
		int verb, int s, int nr_omp_threads,
		bool debug, bool diag, std::string outPath,
		ReferenceMap *reference, ObservationModel *obsModel)
{
	this->verb = verb;
	this->s = s;
	sh = s/2 + 1;
	this->nr_omp_threads = nr_omp_threads;

	this->debug = debug;
	this->diag = diag;
	this->outPath = outPath;

	this->reference = reference;
	this->obsModel = obsModel;

	angpix = obsModel->getPixelSize(0);

	double kmin_px = obsModel->angToPix(kmin, s, 0);
	freqWeight = reference->getHollowWeight(kmin_px);

	ready = true;
}

void BFactorRefiner::processMicrograph(
		long g, MetaDataTable& mdt,
		const std::vector<Image<Complex>>& obs,
		const std::vector<Image<Complex>>& pred)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: BFactorRefiner::processMicrograph: BFactorRefiner not initialized.");
	}

	long pc = obs.size();

	std::stringstream stsg;
	stsg << g;

	std::vector<std::vector<std::pair<int,d2Vector>>> valsPerPart(nr_omp_threads);
	
	const double as = s * angpix;
	const double min_B_px = min_B / (as*as);
	const double max_B_px = max_B / (as*as);
			
	// Parallel loop over all particles in this micrograph
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (long p = 0; p < pc; p++)
	{
		std::stringstream stsp;
		stsp << p;

		CTF ctf;
		ctf.readByGroup(mdt, obsModel, p);
		Image<RFLOAT> ctfImg(sh,s);
		ctf.getFftwImage(ctfImg(), s, s, angpix, false, false, false, false);
		
		Image<Complex> predCTF;
		
		ImageOp::multiply(ctfImg, pred[p], predCTF);
		
		d2Vector sigmaK = BFactorRefiner::findBKRec2D(
			obs[p], predCTF, freqWeight, min_B_px, max_B_px, min_scale, 20, 5, 0.1);
		
		
		int threadnum = omp_get_thread_num();
		
		valsPerPart[threadnum].push_back(std::make_pair(p, sigmaK));
	}
	
	for (int t = 0; t < nr_omp_threads; t++)
	{
		for (int i = 0; i < valsPerPart[t].size(); i++)
		{
			int p = valsPerPart[t][i].first;
			d2Vector BK =  valsPerPart[t][i].second;
			
			if (debug)
			{
				std::cout << p << ": " << as*as*BK[0] << " \t " << BK[1] << "\n";
			}
			
			mdt.setValue(EMDL_CTF_BFACTOR, as*as*BK[0], p);
			mdt.setValue(EMDL_CTF_SCALEFACTOR, BK[1], p);
		}
	}

	// Output a diagnostic Postscript file
	//writeEPS(mdt);

	// Now write out STAR file with optimised values for this micrograph
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	mdt.write(outRoot + "_bfactor_fit.star");
}

std::vector<MetaDataTable> BFactorRefiner::merge(const std::vector<MetaDataTable>& mdts)
{
	int gc = mdts.size();
	int barstep;

	if (verb > 0)
	{
		std::cout << " + Combining data for all micrographs " << std::endl;
		init_progress_bar(gc);
		barstep = 1;
	}

	std::vector<MetaDataTable> mdtOut;
	std::vector<FileName> fn_eps;

	for (long g = 0; g < gc; g++)
	{
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);

		// Read in STAR file with B-factor fit data
		MetaDataTable mdt;
		mdt.read(outRoot+"_bfactor_fit.star");

		mdtOut.push_back(mdt);

		if (exists(outRoot+"_ctf-refine_fit.eps"))
		{
			fn_eps.push_back(outRoot+"_ctf-refine_fit.eps");
		}

		if (verb > 0)
		{
			progress_bar(g);
		}
	}

	if (verb > 0)
	{
		progress_bar(gc);
	}

	if (fn_eps.size() > 0)
	{
		joinMultipleEPSIntoSinglePDF(outPath + "logfile.pdf", fn_eps);
	}
	
	return mdtOut;
}

void BFactorRefiner::writeEPS(const MetaDataTable& mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: BFactorRefiner::writeEPS: BFactorRefiner not initialized.");
	}

	/*std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);

	FileName fn_eps = outRoot + "_defocus_fit.eps";

	CPlot2D plot2D(fn_eps);
	plot2D.SetXAxisSize(600);
	plot2D.SetYAxisSize(600);
	plot2D.SetDrawLegend(false);
	plot2D.SetFlipY(true);

	RFLOAT min_defocus = 99.e10;
	RFLOAT max_defocus = -99.e10;

	const int pc = mdt.numberOfObjects();

	for (int p = 0; p < pc; p++)
	{
		RFLOAT defU, defV;
		mdt.getValue(EMDL_CTF_DEFOCUSU, defU, p);
		mdt.getValue(EMDL_CTF_DEFOCUSV, defV, p);
		defU = (defU + defV) / 2.;

		min_defocus = XMIPP_MIN(min_defocus, defU);
		max_defocus = XMIPP_MAX(max_defocus, defU);
	}

	for (int p = 0; p < pc; p++)
	{
		RFLOAT defU, defV;
		RFLOAT xcoor, ycoor;

		mdt.getValue(EMDL_IMAGE_COORD_X, xcoor, p);
		mdt.getValue(EMDL_IMAGE_COORD_Y, ycoor, p);
		mdt.getValue(EMDL_CTF_DEFOCUSU, defU, p);
		mdt.getValue(EMDL_CTF_DEFOCUSV, defV, p);

		defU = (defU + defV) / 2.;

		RFLOAT val  = (defU - min_defocus) / (max_defocus - min_defocus);

		CDataSet dataSet;

		dataSet.SetDrawMarker(true);
		dataSet.SetDrawLine(false);
		dataSet.SetMarkerSize(10);
		dataSet.SetDatasetColor(val, 0.8 * val, 1.0 - val);

		CDataPoint point(xcoor, ycoor);

		dataSet.AddDataPoint(point);

		plot2D.AddDataSet(dataSet);
	}

	char title[256];
	snprintf(title, 255, "Defocus range from blue to orange: %.0f A", max_defocus - min_defocus);
	plot2D.SetXAxisTitle(title);

	plot2D.OutputPostScriptPlot(fn_eps);*/
}

bool BFactorRefiner::isFinished(const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: BFactorRefiner::isFinished: BFactorRefiner not initialized.");
	}

	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);

	return exists(outRoot + "_bfactor_fit.star");
}

d2Vector BFactorRefiner::findSigmaKRec1D(
		const std::vector<Complex> &obs, 
		const std::vector<Complex> &pred, 
		double sig0, double sig1, 
		int steps, int depth, double rangeItFract)
{
	
}

d2Vector BFactorRefiner::findBKRec2D(
		const Image<Complex> &obs, 
		const Image<Complex> &pred, 
		const Image<RFLOAT> &weight, 
		double B0, double B1, double min_scale,
		int steps, int depth, double rangeItFract)
{
	double minErr = std::numeric_limits<double>::max();
    double bestB = B0;
    double bestA = 1.0;

    const int s = obs.data.xdim;
	const int sh = s/2 + 1;
	
	std::vector<double> sigVals(sh);

    for (int st = 0; st < steps; st++)
    {
        const double B = B0 + st*(B1 - B0)/(steps-1);
		
		for (int r = 0; r < sh; r++)
		{
			sigVals[r] = exp(-B * r * r / 4.0);
		}

        // find optimal scale-factor for hypothetical B-factor
        double num = 0.0, denom = 0.0;

		for (long y = 0; y < s;  y++)
		for (long x = 0; x < sh; x++)
		{
			const int xx = x;
			const int yy = y < sh? y : y - s;
			const int r = (int) (sqrt(xx*xx + yy*yy) + 0.5);
			
			if (r >= sh) continue;
			
			Complex vx = DIRECT_A2D_ELEM(pred.data, y, x);
			const Complex vy = DIRECT_A2D_ELEM(obs.data, y, x);
			const double vw = DIRECT_A2D_ELEM(weight.data, y, x);
			const double vb = sigVals[r];
			
			num   += vw * vb * (vx.real * vy.real + vx.imag * vy.imag);
			denom += vw * vb * vb * (vx.real * vx.real + vx.imag * vx.imag);
		}
		
        const double eps = 1e-20;
        double a = denom > eps? num / denom : num / eps;
		
		if (a < min_scale) a = min_scale;

        double sum = 0.0;

		for (long y = 0; y < s;  y++)
		for (long x = 0; x < sh; x++)
		{
			const int xx = x;
			const int yy = y < sh? y : y - s;
			const int r = (int) (sqrt(xx*xx + yy*yy) + 0.5);
			
			if (r >= sh) continue;
			
			Complex vx = DIRECT_A2D_ELEM(pred.data, y, x);
			const Complex vy = DIRECT_A2D_ELEM(obs.data, y, x);
			const double vw = DIRECT_A2D_ELEM(weight.data, y, x);
			const double vb = sigVals[r];
			
			sum += vw * (vy - a * vb * vx).norm();
		}

        if (sum < minErr)
        {
            minErr = sum;
            bestB = B;
            bestA = a;
        }
    }

    if (depth > 0)
    {
        const double hrange = 0.5 * (B1 - B0);
        double Bnext0 = bestB - rangeItFract*hrange;
        double Bnext1 = bestB + rangeItFract*hrange;

        return findBKRec2D(
			obs, pred, weight, Bnext0, Bnext1, min_scale, 
			steps, depth - 1, rangeItFract);
    }

    return d2Vector(bestB, bestA);
}
