#include <src/jaz/math/fft.h>
#include <src/error.h>
#include <src/image.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/gravis/tImage.h>
#include <src/jaz/image/color_helper.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/interpolation.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;
	
	std::string inPath, outPath, inRef, inMask, inBackgd, inClutter;
	bool has_y_bounds, discr_params;
	double y_min_cl, y_max_cl, ref_scale;
	int amax(72), iters(100), nr_omp_threads(1);
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		inPath = parser.getOption("--i", "Input directory: loads */ref.png, "
													 "*/mask.png and */background.png");		
		inRef = inPath+"/ref.png";
		inMask = inPath+"/mask.png";
		inBackgd = inPath+"/background.png";
		
		std::string refTmp = parser.getOption("--ref", "Reference image (overrides --i)", "");
		std::string maskTmp = parser.getOption("--mask", "Reference mask (overrides --i)", "");
		std::string bgTmp = parser.getOption("--bg", "Background image (overrides --i)", "");
		
		if (refTmp != "") inRef = refTmp;
		if (maskTmp != "") inMask = maskTmp;
		if (bgTmp != "") inBackgd = bgTmp;
		
		iters = textToInteger(parser.getOption("--pc", "Number of particles", "100"));
		amax = textToInteger(parser.getOption("--amax", "Number of angular samples", "360"));
		ref_scale = textToDouble(parser.getOption("--rs", "Reference value scale", "1"));
		
		inClutter = parser.getOption("--cl", "Clutter image: L2 dist. is divided "
		                                     "by the power spectrum of this image", "");
		
		y_min_cl = textToDouble(parser.getOption("--ymin", "Min. y coordinate for particle", "-1"));
		y_max_cl = textToDouble(parser.getOption("--ymax", "Max. y coordinate for particle", "-1"));
		
		has_y_bounds = y_min_cl > 0 && y_max_cl > 0;
		
		discr_params = parser.checkOption("--di", "Discretize unknown positions and angles");
		
		nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
		
		outPath = parser.getOption("--o", "Output path");
				
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	if (outPath[outPath.size() - 1] != '/')
	{
		outPath = outPath + "/";
	}
	
	int res = system(("mkdir -p " + outPath).c_str());
	
	BufferedImage<RFLOAT> ref0, mask, background0;	
	
	ref0.read(inRef);
	ref0 = -ref_scale * ref0;
	
	mask.read(inMask);
		
	background0.read(inBackgd);	
	
	const int wbg = background0.xdim;
	const int hbg = background0.ydim;	
	
	const int wr = ref0.xdim;
	const int hr = ref0.ydim;
	
	const int wb = hbg;
	const int hb = hbg;
	
	
	// center reference value
	
	double avg_rv = 0.0, antimask_sum = 0.0;
	
	for (int y = 0; y < hr; y++)
	for (int x = 0; x < wr; x++)
	{
		avg_rv += (1 - mask(x,y)) * ref0(x,y);
		antimask_sum += (1 - mask(x,y));
	}
	
	avg_rv /= antimask_sum;
	
	ref0 = ref0 - avg_rv;
	
	
	// filter images by spectral power of clutter
	
	BufferedImage<RFLOAT> clutterSpectrum;
	
	if (inClutter != "")
	{	
		BufferedImage<RFLOAT> clutter;
		
		clutter.read(inClutter);
		
		clutterSpectrum = PowerSpectrum::periodogramAverage2D(clutter, wb, hb, 2.0);
		
		Centering::fftwHalfToFftwFull(clutterSpectrum).write(outPath+"powSpectrum.vtk");
	}
	else
	{
		clutterSpectrum = PowerSpectrum::periodogramAverage2D(background0, wb, hb, 2.0);
		
		Centering::fftwHalfToHumanFull(clutterSpectrum).write(outPath+"powSpectrum_2D.vtk");
		
		std::vector<RFLOAT> ra = RadialAvg::fftwHalf_2D_NN(clutterSpectrum);
		
		std::ofstream raof(outPath+"powSpectrum_1D.dat");
		
		for (int i = 1; i < ra.size(); i++)
		{
			raof << i << " " << ra[i] << "\n";
		}
		
		raof.close();
		
		clutterSpectrum = RadialAvg::toFftwHalf_2D_lin(ra, clutterSpectrum.xdim, clutterSpectrum.ydim);
		
		Centering::fftwHalfToHumanFull(clutterSpectrum).write(outPath+"powSpectrum_1D.vtk");
	}
	
	
	BufferedImage<RFLOAT> logClutterPower = Centering::fftwHalfToHumanFull(clutterSpectrum);
	
	for (int y = 0; y < hb; y++)
	for (int x = 0; x < wb; x++)
	{
		if (logClutterPower(x,y) > 0.0)
		{
			logClutterPower(x,y) = log(logClutterPower(x,y));
		}
	}
	
	ColorHelper::writeSignedToPNG(logClutterPower, outPath+"logPowSpectrum.png");	
	
	Normalization::toUnitInterval(background0).write(outPath+"raw_background.png");
	Centering::fftwFullToHumanFull(Normalization::toUnitInterval(ref0)).write(outPath+"raw_ref.png");
	
	BufferedImage<RFLOAT> backgroundFilt = PowerSpectrum::divide(background0, clutterSpectrum, 1e-12);
	BufferedImage<RFLOAT> refFilt = PowerSpectrum::divide(ref0, clutterSpectrum, 1e-12);
	
	
	Normalization::toUnitInterval(backgroundFilt).write(outPath+"filtered_background.png");	
	Centering::fftwFullToHumanFull(Normalization::toUnitInterval(refFilt)).write(outPath+"filtered_ref.png");
	
	
	std::vector<std::string> modeNames{"raw", "filt"};
	std::vector<std::string> metricNames{"CC", "NCC", "wNCC", "wL2"};
	
	
	const double rr = sqrt(wr*wr + hr*hr)/2.0;
	
	const double x_min = rr;
	const double x_max = wbg - rr;
	
	const double y_min = has_y_bounds? y_min_cl : rr;
	const double y_max = has_y_bounds? y_max_cl : hbg - rr;
	
	const int angErrBins = 20, posErrBins = 20, maxPosErr = 200;
	
	// mode      metric      amount
	std::vector<std::vector<std::vector<int>>> angErrHist(2);
	std::vector<std::vector<std::vector<int>>> posErrHist(2);
	
	for (int mode = 0; mode < 2; mode++)
	{
		angErrHist[mode] = std::vector<std::vector<int>>(4);
		posErrHist[mode] = std::vector<std::vector<int>>(4);
		
		for (int met = 0; met < 4; met++)
		{
			angErrHist[mode][met] = std::vector<int>(angErrBins, 0);
			posErrHist[mode][met] = std::vector<int>(posErrBins, 0);
		}
	}
	
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (int it = 0; it < iters; it++)
	{
		double phi_true = 2.0 * PI * rand() / (double)RAND_MAX;
		double x_true = x_min + (x_max - x_min) * rand() / (double)RAND_MAX;
		double y_true = y_min + (y_max - y_min) * rand() / (double)RAND_MAX;

		if (discr_params)
		{
			x_true = (int)(x_true + 0.5);
			y_true = (int)(y_true + 0.5);
			
			const int phiInd = (int) (amax * phi_true / (2.0 * PI) + 0.5);
			phi_true = 2.0 * PI * phiInd / (double)amax;
		}
		
		const d2Matrix R(cos(phi_true),-sin(phi_true),
						 sin(phi_true), cos(phi_true));
		
		BufferedImage<RFLOAT> data0 = background0;
		
		const d2Vector centref(wr / 2, hr / 2);
		
		for (int y = 0; y < hbg; y++)
		for (int x = 0; x < wbg; x++)
		{
			const d2Vector pc(x - x_true, y - y_true);
			const d2Vector pr = centref + R * pc;
			
			if (pr.x < 0 || pr.x >= wr || pr.y < 0 || pr.y >= hr)
			{
				continue;
			}
			
			const RFLOAT v = Interpolation::linearXY_clip(ref0, pr.x, pr.y);
			const RFLOAT m = Interpolation::linearXY_clip(mask, pr.x, pr.y);
			
			data0(x,y) += m*v;
		}
		
		for (int mode = 0; mode < 2; mode++)
		{
			BufferedImage<RFLOAT> data;
			
			if (mode == 1)
			{
				data = PowerSpectrum::divide(data0, clutterSpectrum, 1e-12);
			}
			else
			{
				data = data0;
			}
			
			if (it == 0)
			{
				Normalization::toUnitInterval(data).write(outPath+"data_"+modeNames[mode]+".png");
			}
			
			data = Centering::humanFullToFftwFull(data);

			std::vector<double> bestVal(4, -std::numeric_limits<double>::max()), 
					bestAng(4, -1.0), bestX(4, 0.0), bestY(4, 0.0);
					
			for (int a = 0; a < amax; a++)
			{
				const double phi_hyp = 2.0 * PI * a / (double)amax;
				
				const d2Matrix R2(cos(phi_hyp),-sin(phi_hyp),
								  sin(phi_hyp), cos(phi_hyp));
				
				BufferedImage<RFLOAT> refRot(wr,hr), maskRot(wr, hr);
				
				for (int y = 0; y < hr; y++)
				for (int x = 0; x < wr; x++)
				{
					const d2Vector pc(x - centref.x, y - centref.y);
					const d2Vector pr = centref + R2 * pc;
					
					if (pr.x < 0 || pr.x >= wr || pr.y < 0 || pr.y >= hr)
					{
						continue;
					}
					
					refRot(x,y) = Interpolation::linearXY_clip(ref0, pr.x, pr.y);
					maskRot(x,y) = Interpolation::linearXY_clip(mask, pr.x, pr.y);
				}
				
				refRot = Centering::humanFullToFftwFull(refRot);
				maskRot = Centering::humanFullToFftwFull(maskRot);
				
				if (mode == 1)
				{
					refRot = PowerSpectrum::divide(refRot, clutterSpectrum, 1e-12);
				}
						
				std::vector<BufferedImage<RFLOAT>> metrics(4);
				
				metrics[0] = Centering::fftwFullToHumanFull(Similarity::CC_2D(refRot, data, true));
				metrics[1] = Centering::fftwFullToHumanFull(Similarity::NCC_2D(refRot, data, true));
				metrics[2] = Centering::fftwFullToHumanFull(Similarity::weightedNCC_2D(refRot, maskRot, data, 1e-4));
				metrics[3] = -Centering::fftwFullToHumanFull(Similarity::weightedL2_2D(refRot, maskRot, data));
				
				/*{
					std::stringstream sts;
					sts << phi_hyp;
					
					metrics[3].write(outPath+"mL2_"+modeNames[mode]+"_phi-"+sts.str()+".vtk");
				}*/
				
				//std::cout << phi_hyp << ": ";
							 
				for (int met = 0; met < 4; met++)
				{
					d3Vector optPosVal = discr_params?
						Interpolation::discreteMaxXYV(metrics[met]) :
						Interpolation::quadraticMaxXYV(metrics[met]);
					
					//std::cout << optPosVal.z << "  ";
								 
					if (optPosVal.z > bestVal[met])
					{
						bestVal[met] = optPosVal.z;
						bestAng[met] = phi_hyp;
						bestX[met] = optPosVal.x;
						bestY[met] = optPosVal.y;
					}
					
					/*if (a == 0) 
					{
						if (met < 3)
						{
							ColorHelper::writeSignedToPNG(
								metrics[met], 
								outPath + "sim_" + metricNames[met] + "_" + modeNames[mode] + ".png");
						}
						else
						{
							ColorHelper::writeSignedToPNG(
								Normalization::toUnitInterval(metrics[met]), 
								outPath + "sim_" + metricNames[met] + "_" + modeNames[mode] + ".png");
						}
					}*/
				}
				
				//std::cout << std::endl;
			}
			
			#pragma omp critical	
			{
				std::cout << "particle " << it << "\n";
				std::cout << "phi_true = " << phi_true << std::endl;
				std::cout << "x/y_true = " << x_true << ", " << y_true << std::endl;
				
				std::cout << "ang. error (" << modeNames[mode] << "): ";
			
				for (int met = 0; met < 4; met++)
				{
					double angErr = bestAng[met] - phi_true;
					double angErr_abs = std::abs(angErr);
					
					if (angErr_abs > PI) angErr_abs = 2.0*PI - angErr_abs;
					
					int angErrInd = (int)(angErrBins * angErr_abs / PI + 0.5);
					
					if (angErrInd >= angErrBins)
					{
						angErrInd = angErrBins - 1;
					}
					
					angErrHist[mode][met][angErrInd]++; 
						
					std::cout << angErr << " (";
					
					d2Vector posErr(bestX[met] - x_true, bestY[met] - y_true);
					double posErr_abs = posErr.length();
					
					std::cout << posErr_abs << ")  ";
					
					int posErrInd = (int)(posErrBins * posErr_abs / maxPosErr);
					
					if (posErrInd >= posErrBins)
					{
						posErrInd = posErrBins - 1;
					}
					
					posErrHist[mode][met][posErrInd]++;
				}
				
				std::cout << std::endl;
			}
			
			
			
			/*Image<RFLOAT> CC = Similarity::CC_2D(refFilt, backgroundFilt, true);
			ColorHelper::writeSignedToPNG(Centering::fftwFullToHumanFull(CC), outPath+"sim_CC.png");
			
			Image<RFLOAT> NCC = Similarity::NCC_2D(refFilt, backgroundFilt, true);
			ColorHelper::writeSignedToPNG(Centering::fftwFullToHumanFull(NCC), outPath+"sim_NCC.png");	
			
			Image<RFLOAT> wNCC = Similarity::weightedNCC_2D(refFilt, mask, backgroundFilt, 1e-4);
			ColorHelper::writeSignedToPNG(Centering::fftwFullToHumanFull(wNCC), outPath+"sim_wNCC.png");
			
			Image<RFLOAT> wL2 = Similarity::weightedL2_2D(refFilt, mask, backgroundFilt);
			ColorHelper::writeSignedToPNG(Centering::fftwFullToHumanFull(wL2), outPath+"sim_wL2.png");*/
		}
	}
	
	for (int mode = 0; mode < 2; mode++)
	{
		for (int met = 0; met < 4; met++)
		{
			std::ofstream posHistFile(outPath + "posErrHist_" 
				+ metricNames[met] + "_" + modeNames[mode] + ".dat");
			
			const double posHistScale = maxPosErr / (double) posErrBins;
					
			for (int i = 0; i < posErrBins; i++)
			{
				posHistFile << i * posHistScale << " " << posErrHist[mode][met][i] << "\n";
				posHistFile << (i+1) * posHistScale << " "  << posErrHist[mode][met][i] << "\n";
			}
			
			
			std::ofstream angHistFile(outPath + "angErrHist_" 
				+ metricNames[met] + "_" + modeNames[mode] + ".dat");
			
			const double angHistScale = PI / (double) angErrBins;
					
			for (int i = 0; i < angErrBins; i++)
			{
				angHistFile << i * angHistScale << " "  << angErrHist[mode][met][i] << "\n";
				angHistFile << (i+1) * angHistScale << " "  << angErrHist[mode][met][i] << "\n";
			}
		}
	}
	
	return 0;
}
