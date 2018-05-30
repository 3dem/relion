#include "bfactor_estimator.h"
#include "motion_refiner.h"
#include "motion_helper.h"

#include <src/jaz/micrograph_handler.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/image_log.h>
#include <src/filename.h>
#include <src/backprojector.h>

using namespace gravis;

BFactorEstimator::BFactorEstimator()
{}

void BFactorEstimator::read(IOParser& parser, int argc, char* argv[])
{
    parser.addSection("Old B-factor estimation options");

    doAnything = parser.checkOption("--old_bfactors", "Estimate classical B/C-factors.");
	f_max = textToInteger(parser.getOption("--old_bfac_f_max", "Number of frames to consider", "-1"));
	f_batch = textToInteger(parser.getOption("--old_bfac_f_batch", "Number of frames per batch", "10"));
			
    bfac_diag = parser.checkOption("--diag_old_bfactor", "Write out old B/C-factor diagnostic data");
	
	fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
	nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
	helical_rise = textToDouble(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
	helical_twist = textToDouble(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));
	
	fit_minres = textToDouble(parser.getOption("--autob_lowres", "Lowest resolution (in A) to include in fitting of the B-factor", "20."));
	perframe_highres = textToDouble(parser.getOption("--perframe_highres", "Highest resolution (in A) to include in the per-frame reconstructions", "-1."));

}

void BFactorEstimator::init(
    int verb, int s, int fc, 
    int nr_omp_threads, std::string outPath, bool debug,
    ObservationModel* obsModel,
    MicrographHandler* micrographHandler,
	ReferenceMap* reference)
{
    this->verb = verb;
    this->s = s;
    this->sh = s/2 + 1;
    this->fc = fc;
    this->nr_omp_threads = nr_omp_threads;
    this->outPath = outPath;
    this->debug = debug;
    this->obsModel = obsModel;
    this->micrographHandler = micrographHandler;
    this->angpix = obsModel->angpix;
	this->reference = reference;
}

void BFactorEstimator::process(const std::vector<MetaDataTable>& mdts)
{	
	if (f_max < 0) f_max = fc-1;
	
	const int bc = (int) ceil((f_max + 1) / (double)f_batch);			
	const int gc = mdts.size();
	
	for (int b = 0; b < bc; b++)
	{
		std::cout << " + Batch " << (b+1) << " of " << bc << std::endl;
		const int f0 = b * f_batch;
		int f1 = (b+1) * f_batch - 1;
		
		if (f1 > f_max) f1 = f_max;
		
		const int fcb = f1 - f0 + 1;
	
		int barstep;
		int my_nr_micrographs = gc;
	
		if (verb > 0)
		{
			std::cout << " +    Estimating classical B-factors ... " << std::endl;
		}
	
		std::vector<ParFourierTransformer> fts(nr_omp_threads);
	
		long nr_done = 0;
		
		if (verb > 0)
		{
			std::cout << " +       Allocating " << (2*fcb) << " backprojectors ... " << std::endl;
		}
		
		std::vector<BackProjector> backprojectors(2*fcb);
		
		for (int f = 0; f < fcb; f++)
		{
			for (int h = 0; h < 2; h++)
			{
				backprojectors[2*f+h] = BackProjector(s, 3, fn_sym, TRILINEAR, 1.0);
				backprojectors[2*f+h].initZeros(s);
			}
		}
			
		if (verb > 0)
		{
			std::cout << " +       Backprojecting frames ... " << std::endl;
			init_progress_bar(my_nr_micrographs);
			barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
		}
		
		Image<RFLOAT> ctfImg(sh,s), ctfImgSq(sh,s);
	
		for (long g = 0; g < gc; g++)
		{
			const int pc = mdts[g].numberOfObjects();
	
			std::vector<std::vector<Image<Complex>>> movie;
			movie = micrographHandler->loadMovie(mdts[g], s, angpix, fts);
	
			FileName fn_root = MotionRefiner::getOutputFileNameRoot(outPath, mdts[g]);
			std::vector<std::vector<d2Vector>> shift;
			shift = MotionHelper::readTracks(fn_root+"_tracks.star");
	
			for (int p = 0; p < pc; p++)
			{
				int subset;
				mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, subset);
				subset -= 1;
				
				CTF ctf;
				ctf.read(mdts[g], mdts[g], p);
				ctf.getFftwImage(ctfImg(), s, s, angpix, false, false, false, true);
				
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(ctfImgSq())
				{
					DIRECT_MULTIDIM_ELEM(ctfImgSq(), n) = 
							DIRECT_MULTIDIM_ELEM(ctfImg(), n) * DIRECT_MULTIDIM_ELEM(ctfImg(), n);
				}
				
				RFLOAT rot, tilt, psi;
				mdts[g].getValue(EMDL_ORIENT_ROT, rot, p);
				mdts[g].getValue(EMDL_ORIENT_TILT, tilt, p);
				mdts[g].getValue(EMDL_ORIENT_PSI, psi, p);
				
				Matrix2D<RFLOAT> A3D;
				Euler_angles2matrix(rot, tilt, psi, A3D);
				
				Matrix1D<RFLOAT> trans(2);
				trans.initZeros();
				mdts[g].getValue(EMDL_ORIENT_ORIGIN_X, XX(trans), p);
				mdts[g].getValue(EMDL_ORIENT_ORIGIN_Y, YY(trans), p);
				
				
				#pragma omp parallel for num_threads(nr_omp_threads)
				for (int f = 0; f < fcb; f++)
				{
					const int ff = f0 + f;
					
					Image<Complex> obs(sh,s);
					
					shiftImageInFourierTransform(
						movie[p][ff](), obs(), s, 
						XX(trans) - shift[p][ff].x - s/2, 
						YY(trans) - shift[p][ff].y - s/2);
									
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(obs())
					{
						DIRECT_MULTIDIM_ELEM(obs(), n) *= DIRECT_MULTIDIM_ELEM(ctfImg(), n);
					}
	
					backprojectors[2*f + subset].set2DFourierTransform(
								obs(), A3D, IS_NOT_INV, &ctfImgSq());
				}
			}
	
			nr_done++;
	
			if (verb > 0 && nr_done % barstep == 0)
			{
				progress_bar(nr_done);
			}
		}

		if (verb > 0)
		{
			progress_bar(my_nr_micrographs);
		}
		
		if (verb > 0)
		{
			std::cout << " +       Performing reconstruction ... " << std::endl;
		}
		
		//#pragma omp parallel for num_threads(nr_omp_threads)
		for (int f = 0; f < fcb; f++)
		{
			const int ff = f0 + f;
			
			for (int h = 0; h < 2; h++)
			{
				backprojectors[2*f + h].symmetrise(nr_helical_asu, helical_twist, helical_rise / angpix);
			
				MultidimArray<RFLOAT> dummy;
				Image<RFLOAT> vol;
				backprojectors[2*f + h].reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy, dummy);
			
				std::string fnOut = outPath + "frame_" + integerToString(ff+1, 3, '0')
					+ "_half_" + integerToString(h+1) + ".mrc";
				
				vol.write(fnOut);
			}
		}
	}
	
	if (verb > 0)
    {
        std::cout << " +    Computing B-factors ... " << std::endl;
	}
	
	std::vector<Image<RFLOAT>> avgs(2);
	
	for (int h = 0; h < 2; h++)
	{
		avgs[h] = Image<RFLOAT>(s,s,s);
		
		for (int f = 0; f <= f_max; f++)
		{
			std::string fnIn = outPath + "frame_" + integerToString(f+1, 3, '0')
				+ "_half_" + integerToString(h+1) + ".mrc";
			
			Image<RFLOAT> vol;
			vol.read(fnIn);
			avgs[h]() += vol();
		}
		
		std::string fnAvgOut = outPath + "avg_half_" + integerToString(h+1) + ".mrc";		
		avgs[h].write(fnAvgOut);
	}
	
	MultidimArray<RFLOAT> avgFSC;
	
	if (reference->hasMask) 
	{
		avgFSC = maskedFSC(avgs[0], avgs[1], reference->mask);
	}
	else 
	{
		getFSC(avgs[0](), avgs[1](), avgFSC);
	}
	
	MultidimArray<RFLOAT> perframe_bfactors(3*(f_max+1));
			
	std::vector<Image<RFLOAT>> vols(2);
	
	for (int f = 0; f <= f_max; f++)
	{
		for (int h = 0; h < 2; h++)
		{
			std::string fnIn = outPath + "frame_" + integerToString(f+1, 3, '0')
				+ "_half_" + integerToString(h+1) + ".mrc";
			
			vols[h].read(fnIn);
		}
		
		MultidimArray<RFLOAT> fsc;
		
		if (reference->hasMask) 
		{
			fsc = maskedFSC(vols[0], vols[1], reference->mask);
		}
		else 
		{
			getFSC(vols[0](), vols[1](), fsc);
		}
		
		double bfactor, offset, corr_coeff;
		
		calculateBfactorSingleFrameReconstruction(
			f, fsc, avgFSC, bfactor, offset, corr_coeff);
		
		DIRECT_A1D_ELEM(perframe_bfactors, f * 3 + 0) = bfactor;
		DIRECT_A1D_ELEM(perframe_bfactors, f * 3 + 1) = offset;
		DIRECT_A1D_ELEM(perframe_bfactors, f * 3 + 2) = corr_coeff;
	}
	
	writeStarFileBfactors(perframe_bfactors);
}

bool BFactorEstimator::doingAnything()
{
	return doAnything;
}

MultidimArray<RFLOAT> BFactorEstimator::maskedFSC( 
		Image<RFLOAT>& I1, 
		Image<RFLOAT>& I2,
		const Image<RFLOAT>& Imask)
{
	MultidimArray<RFLOAT> out;
	
	I1() *= Imask();
	I2() *= Imask();
	
	getFSC(I1(), I2(), out);
	
	return out;
}

bool BFactorEstimator::calculateBfactorSingleFrameReconstruction(
		int frame,
		const MultidimArray<RFLOAT>& fsc_frame,
		const MultidimArray<RFLOAT>& fsc_average,
		double& bfactor, double& offset, double& corr_coeff)
{
	// Now use relative ratio of signal amplitudes w.r.t. the average of all single-frame reconstructions
	// SSNR= tau^2/sigma^2 = FSC / (1 - FSC)
	// tau_frame / tau_avg = tau_f / tau_a = sqrt (FSC_f / (1 - FSC_f)) / sqrt (FSC_a / (1 - FSC_a))
	// = sqrt( {FSC_f / (1 - FSC_f)} * {(1 - FSC_a) / FSC_a} )
	// = sqrt( (FSC_f - FSC_f * FSC_a) / (FSC_a - FSC_f * FSC_a)  )
	// Then make a Guinier plot of that: ln(tau) versus 1/d^2 and
	// fit a line through all points where FSC_f < 1 && FSC_f > 0.143
	// Store the bfactor (4*slope) AND the offset of that line

	MetaDataTable MDout;
	MDout.setName("relative_guinier");

	fit_point2D onepoint;
	std::vector<fit_point2D> guinier;
	
	for (int i = 1; i < XSIZE(fsc_frame); i++) // ignore origin
	{
		RFLOAT res = (RFLOAT)i / (s * angpix); // resolution in 1/A
		RFLOAT resang = 1. / res;
		RFLOAT res2 = res*res;

		RFLOAT fsc_f = DIRECT_A1D_ELEM(fsc_frame, i);
		RFLOAT fsc_a = DIRECT_A1D_ELEM(fsc_average, i);

		if (resang < fit_minres && resang > perframe_highres 
				&& fsc_f < 1 && fsc_a < 1 && fsc_f > 0.143 && fsc_a > 0.143)
		{
			// ln(tau_f / tau_a)
			// I could have calculated the same using: ln(tau_f / tau_a)   = ln(tau_f) - ln(tau_a)
			// where tau_f = sqrt (FSC_f / (1 - FSC_f)) and tau_a = sqrt (FSC_a / (1 - FSC_a))
			// This is numerically identical
			RFLOAT logreltau = log( sqrt( (fsc_f - fsc_f * fsc_a) / (fsc_a - fsc_f * fsc_a) ) );

			onepoint.x = res2;
			onepoint.y = logreltau;
			onepoint.w = 1.;
			guinier.push_back(onepoint);

			MDout.addObject();
			MDout.setValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, res2);
			MDout.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_IN, logreltau);
		}
	}
	
	MDout.write(outPath + "frame_" + integerToString(frame+1, 3, '0') + "_guinier.star");

	// Check if any points were included in the Guinier plot
	if (guinier.size() < 3)
	{
		std::cerr << " WARNING: insufficient number of points in the Guinier "
		          << "plot of movie frame: " << frame << std::endl;
		
		std::cerr << " Consider lowering the lowres-limit, or average over "
		          << "multiple frames in the B-factor estimation." << std::endl;
		
		bfactor = 0.0;
		offset = 0.0;
		corr_coeff = 0.0;
		
		return false;
	}

	// Now do the fit
	fitStraightLine(guinier, bfactor, offset, corr_coeff);
	// this is the B-factor relative to the average from all single-frame reconstructions!
	// in this case: positive values mean BETTER than average, and thus HIGHER WEIGHTS!
	bfactor *= 4.;

	CDataSet dataSet;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDout)
	{
		RFLOAT res2;
		MDout.getValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, res2);
		CDataPoint point(res2, offset + (res2 * bfactor / 4.));
		dataSet.AddDataPoint(point);
	}
	
	dataSet.SetDatasetColor(1., 0., 0.);
	dataSet.SetDatasetTitle("Fitted straight line");
	
	CPlot2D* plot2D = new CPlot2D("Guinier plot frame " + integerToString(frame+1));
	
	plot2D->SetXAxisSize(600);
	plot2D->SetYAxisSize(400);
	plot2D->SetXAxisTitle("resolution (1/A^2)");
	plot2D->SetYAxisTitle("ln(amplitudes)");
	MDout.addToCPlot2D(plot2D, EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, EMDL_POSTPROCESS_GUINIER_VALUE_IN);
	plot2D->AddDataSet(dataSet);
	plot2D->OutputPostScriptPlot(outPath + "frame_"+integerToString(frame+1, 3, '0')+"_guinier.eps");
	
	delete plot2D;
	
	return true;
}

void BFactorEstimator::writeStarFileBfactors(
		MultidimArray<RFLOAT>& perframe_bfactors)
{
	MetaDataTable MDout;
	MDout.setName("perframe_bfactors");

	for (int iframe = 0; iframe <= f_max; iframe++)
	{
		MDout.addObject();
		MDout.setValue(EMDL_IMAGE_FRAME_NR, iframe);
		MDout.setValue(EMDL_POSTPROCESS_BFACTOR, DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0) );
		MDout.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1) );
	}

	MDout.write(outPath+"old_bfactors.star");

	CPlot2D *plot2D=new CPlot2D("Polishing B-factors");
	plot2D->SetXAxisSize(600);
	plot2D->SetYAxisSize(400);
	plot2D->SetDrawLegend(false);
	plot2D->SetXAxisTitle("movie frame");
	plot2D->SetYAxisTitle("B-factor");
	MDout.addToCPlot2D(plot2D, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_BFACTOR);
	plot2D->OutputPostScriptPlot(outPath + "old_bfactors.eps");
	delete plot2D;

	CPlot2D *plot2Db=new CPlot2D("Polishing scale-factors");
	plot2Db->SetXAxisSize(600);
	plot2Db->SetYAxisSize(400);
	plot2Db->SetDrawLegend(false);
	plot2Db->SetXAxisTitle("movie frame");
	plot2Db->SetYAxisTitle("Scale-factor");
	MDout.addToCPlot2D(plot2Db, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT);
	plot2Db->OutputPostScriptPlot(outPath + "old_scalefactors.eps");
	delete plot2Db;

}
