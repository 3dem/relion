#include "tomo_ctf_find.h"
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/filter.h>
#include <src/ctf.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/astigmatism_fit.h>
#include <src/jaz/optics/ctf_equiphase_fit.h>
#include <src/jaz/optics/spectral_ctf_cost.h>
#include <src/jaz/optics/tile_ctf_cost.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/util/zio.h>

using namespace gravis;

#define DAMAGE_EPS 1e-12


TomoCtfFind::TomoCtfFind(
	const BufferedImage<float> &tiltSeries, 
	const std::vector<gravis::d4Matrix> &projections, 
	bool diag,
	double pixelSize, double voltage, double Cs, double Q0, 
	int tileSize, double r0_ang, double r1_ang, int firstFrame)
	:	tiltSeries(tiltSeries), 
		projections(projections), 
		diag(diag),
		pixelSize(pixelSize), 
		voltage(voltage), 
		Cs(Cs), 
		Q0(Q0),
		lambda(12.2643247 / sqrt(voltage * 1e3 * (1. + voltage * 0.978466e-3))),
		tileSize(tileSize),
		firstFrame(firstFrame),
		r0_ang(r0_ang),
		r1_ang(r1_ang)
{	
}

void TomoCtfFind::initialize(int num_threads)
{
	this->num_threads = num_threads;
	
	const int w = tiltSeries.xdim;
	const int h = tiltSeries.ydim;
	
	const int s = tileSize;
	const int sh = s/2 + 1;
	
	fc = tiltSeries.zdim;
	
	r0_pix = tileSize * pixelSize / r0_ang;
	r1_pix = tileSize * pixelSize / r1_ang;
	
	if (r0_pix < 1) r0_pix = 1;
	if (r1_pix >= sh) r1_pix = sh - 1;
	
	r0ref_pix = 2 * r0_pix;
	r1ref_pix = sh - 1;
			
	const double overlap = 2.0;
	
	tw = (int)(overlap * (w-s) / s + 0.5);
	th = (int)(overlap * (h-s) / s + 0.5);
	
	tdx = (w - s) / (double) tw;
	tdy = (h - s) / (double) th;
	
	tiles.resize(sh, s, tw*th*fc);
	tileOffZ.resize(tw,th,fc);
	radAvg.resize(tw,th,fc);	
	radVar.resize(tw,th,fc);
	radAvgNrm.resize(tw,th,fc);
	
	
	if (diag)
	{
		ZIO::makeDir("debug");
	}
	
	
	frameDose = Damage::estimateDoseWeights(
		tiltSeries, pixelSize, tileSize, firstFrame, 20.0, diag);
	
	damageWeight.resize(fc);
	
	for (int f = 0; f < fc; f++)
	{
		damageWeight[f] = Damage::weightVector(frameDose[f], pixelSize, tileSize);
	}
	
	frqWgh.resize(sh);
				
	for (int r = 0; r < sh; r++)
	{
		frqWgh[r] = r / (double)sh;
	}
	
	std::cout << "    extracting tiles..." << std::endl;
	
	extractTiles();
	
	std::cout << "    averaging tiles..." << std::endl;
	
	averageTiles();
	
	std::cout << "    subtracting backgrounds..." << std::endl;
	
	subtractBackground();
	
	std::cout << "    computing radial averages..." << std::endl;
	
	averageRadially();
}

void TomoCtfFind::findInitialCtf(double z0, double z1, double dz)
{
	const int s = tileSize;
	const int sh = s/2 + 1;
	const int zc = (int)((z1 - z0)/dz);
	
	std::vector<double> zvals(zc);
	
	for (int i = 0; i < zc; i++)
	{
		zvals[i] = z0 + i * dz;
	}
	
	std::vector<double> costs(zc, 0.0);
	
	double minCost(std::numeric_limits<double>::max());
	double bestZ(0.0);
	
	std::cout << "    scanning for optimal global defoci between " 
			  << z0 << " Å and " << z1 << " Å..." << std::endl;
			
	for (int zi = 0; zi < zc; zi++)
	{
		const double z = zvals[zi];	
		
		double cost(0.0);
		
		CTF ctf;
		ctf.setValues(z, z, 0.0, voltage, Cs, Q0, 0.0);
		
		for (int r = r0_pix; r <= r1_pix; r++)
		{
			const double ra = r / (double)(s * pixelSize);
			
			const double c = ctf.getCTF(ra, 0.0);
			const double d = radAvgAllNrm[r];
			
			cost -= frqWgh[r] * (c*c - 0.5) * d;
		}
		
		costs[zi] = cost;
		
		if (cost < minCost)
		{
			minCost = cost;
			bestZ = z;
		}
	}
	
	if (diag)
	{
		std::ofstream ofs("debug/cost_over_z_all.dat");
		
		for (int zi = 0; zi < zc; zi++)
		{
			ofs << zvals[zi] << " " << costs[zi] << "\n";	
		}
	}
	
	SpectralCtfCost scf(
		avgTile, avgBackground, pixelSize, voltage, Cs, Q0, 
		num_threads, r0ref_pix, r1ref_pix);
	
	std::vector<double> x0 = scf.getInitialParams(bestZ);
	std::vector<double> x1 = LBFGS::optimize(x0, scf, diag? 2 : 0, 2000);
	
	
	if (diag)
	{
		BufferedImage<float> fitImg = scf.render(x1);
		scf.render(x0).write("debug/global_fit_0.mrc");
		fitImg.write("debug/global_fit_1.mrc");
		
		std::vector<float> raFit = RadialAvg::fftwHalf_2D_lin(fitImg, 2);
		std::vector<float> raData = RadialAvg::fftwHalf_2D_lin(avgTile, 2);
		
		ZIO::writeDat(raFit, "debug/global_fit.dat");
		ZIO::writeDat(raData, "debug/global_data.dat");
		
		avgTile.write("debug/global_fit_tgt.mrc");
	}
	
	globalParams = x1;
}

std::vector<double> TomoCtfFind::findDefoci(double z0, double z1, double dz)
{
	const int zc = (int)((z1 - z0)/dz);
	
	std::vector<double> zvals(zc);
	
	for (int i = 0; i < zc; i++)
	{
		zvals[i] = z0 + i * dz;
	}
	
	std::vector<std::vector<std::vector<double>>> cost(2);
	std::vector<std::vector<double>> bestZcost(2);
	std::vector<std::vector<double>> bestZ_perFrame(2);
	
	std::cout << "    scanning for optimal defoci between " 
			  << z0 << " Å and " << z1 << " Å..." << std::endl;
			
	/*for (int hand = -1; hand <= 1; hand++)
	{
		const int hi = (hand + 1)/2;
		
		cost[hi] = std::vector<std::vector<double>>(fc);
		bestZcost[hi] = std::vector<double>(fc);
		bestZ_perFrame[hi] = std::vector<double>(fc);*/
		
		//#pragma omp parallel for num_threads(num_threads)
		for (int f = 20; f < fc; f++)
		{			
			//cost[hi][f] = std::vector<double>(zc,0.0);
			
			RawImage<float> tilesSlab = tiles.getSlabRef(f * tw*th, tw*th);			
			RawImage<double> zOffsets = tileOffZ.getSliceRef(f);
					
			TileCtfCost tcf(
				tilesSlab, tw, th,
				zOffsets, avgBackground, pixelSize, voltage, Cs, Q0,
				-1.0,
				num_threads, r0ref_pix, r1ref_pix);
			
			std::vector<double> globParamsExp = tcf.getInitialParams(globalParams);
			
			std::vector<double> x1 = LBFGS::optimize(globParamsExp, tcf, diag? 2 : 0, 100);
			
			if (diag)
			{
				std::stringstream sts;
				sts << f;
				
				tcf.renderAlignedSpectrum(x1).write("debug/frame_"+sts.str()+"_fit-central.mrc");
				tcf.renderAvgSpectrum().write("debug/frame_"+sts.str()+"_tgt-central.mrc");
				tcf.renderTable(x1).write("debug/frame_"+sts.str()+"_fit-table.mrc");				
				
				ZIO::writeDat(tcf.radialAverageAligned(x1), 
							  "debug/frame_"+sts.str()+"_radial-avg-ali.dat");
				
				ZIO::writeDat(tcf.radialAverage(x1), 
							  "debug/frame_"+sts.str()+"_radial-avg.dat");
				
				ZIO::writeDat(tcf.radialAverageFit(x1), 
							  "debug/frame_"+sts.str()+"_radial-avg-fit.dat");
						
				/*std::vector<float> raFit = RadialAvg::fftwHalf_2D_lin(fitImg, 2);
				std::vector<float> raData = RadialAvg::fftwHalf_2D_lin(avgTiles[f], 2);
				
				ZIO::writeDat(raFit, "debug/frame_"+sts.str()+"_fit.dat");
				ZIO::writeDat(raData, "debug/frame_"+sts.str()+"_tgt.dat");*/
			}
		}
		
		return std::vector<double>(0);
	//}
	
	std::vector<double> minCostSumByHand(2, 0.0);
	
	for (int hi = 0; hi < 2; hi++)
	{
		for (int f = 0; f < fc; f++)
		{
			minCostSumByHand[hi] += bestZcost[hi][f];
		}
	}
	
	int likelyHand(0);
	
	if (minCostSumByHand[0] < minCostSumByHand[1]) likelyHand = 0;
	else likelyHand = 1;
	
	std::vector<double> out(fc, 0.0);
	
	for (int f = 0; f < fc; f++)
	{
		out[f] = bestZ_perFrame[likelyHand][f];
	}
	
	if (diag)
	{
		std::ofstream ofsH0("debug/CTF_cost_over_z_good-hand.dat");
		std::ofstream ofsH1("debug/CTF_cost_over_z_bad-hand.dat");
		
		for (int f = 0; f < fc; f++)
		{
			std::stringstream sts;
			sts << f;
			
			for (int zi = 0; zi < zc; zi++)
			{
				ofsH0 << (z0 + zi * dz) << " " << cost[likelyHand][f][zi] << "\n";
				ofsH1 << (z0 + zi * dz) << " " << cost[(likelyHand+1)%2][f][zi] << "\n";
			}
			
			ofsH0 << "\n";
			ofsH1 << "\n";
			
			std::ofstream ofs("debug/CTF_cost_over_z_"+sts.str()+".dat");
			
			for (int hand = -1; hand <= 1; hand+=2)
			{
				const int hi = (hand + 1)/2;
				
				for (int zi = 0; zi < zc; zi++)
				{
					ofs << (z0 + zi * dz) << " " << cost[hi][f][zi] << "\n";
				}
				
				ofs << "\n";
			}
			
			std::ofstream ofs2("debug/CTF_1D_fit_"+sts.str()+".dat");
			
			for (int ty = 0; ty < th; ty++)
			for (int tx = 0; tx < tw; tx++)	
			{
				for (int r = r0_pix; r <= r1_pix; r++)
				{
					const double v = radAvgNrm(tx, ty, f)[r];
					ofs2 << r << " " << v << "\n";
				}
				
				ofs2 << "\n";
			}
			
			for (int hi = 0; hi < 2; hi++)
			{			
				const double z = bestZ_perFrame[hi][f];
						
				CTF ctf;
				ctf.setValues(z, z, 0.0, voltage, Cs, Q0, 0.0);
				
				for (int r = r0_pix; r <= r1_pix; r++)
				{
					double ra = r / (double)(tileSize * pixelSize);
					double c = ctf.getCTF(ra, 0.0);
					
					ofs2 << r << " " << c*c << "\n";
				}
				
				ofs2 << "\n";
			}
			
		}
	}
	
	return out;
}

void TomoCtfFind::findAstigmatism()
{
	const int s = tileSize;
	const int sh = s/2 + 1;
	
	const int r0_pix_ast = tileSize * pixelSize / 25.0;
	const int r1_pix_ast = tileSize * pixelSize / 10.0;
	
	
	/*{
		CtfEquiphaseFit cef(
			avgTiles[20], pixelSize, voltage, Cs, Q0,
			0.1, num_threads, r0_pix_ast, r1_pix_ast);
		
		CTF ctf0;
		ctf0.setValues(20000, 22000, 60, voltage, Cs, Q0, 0.0, 1.0, 0.0);
				
		std::vector<double> initial3(3);
		initial3[0] = ctf0.getAxx();
		initial3[1] = ctf0.getAxy();
		initial3[2] = ctf0.getAyy();
		
		std::vector<double> eqFit = NelderMead::optimize(
			initial3, cef, 0.01, 0.000001, 500, 1.0, 2.0, 0.5, 0.5, true);
			
		CTF ctf2 = cef.paramsToCtf(eqFit);
		const int rad = 2 * sh;		
		std::vector<double> accum(rad), weight(rad);
		
		cef.averageAlongIso(ctf2, accum, weight, 2);
		
		Image<double> equiphaseAvg = cef.computeExpansion(ctf2, accum, 2);		
				
		avgTiles[20].write("debug/eq_data.mrc");
		equiphaseAvg.write("debug/eq_fit.mrc");

		std::exit(0);		
	}*/
	
	const bool global_astigmatism = false;
	
	const double Cs_px = (PI / 2.0) * 1e7 * Cs * lambda * lambda * lambda;
			
	AstigmatismFit afGlob(avgTiles, num_threads, r0_pix_ast, r1_pix_ast, Cs_px);
	std::vector<double> initial(2, 0.0);
	
	double alphaGlob, betaGlob;
	
	if (global_astigmatism)
	{
		std::vector<double> globalFit = NelderMead::optimize(
			initial, afGlob, 0.01, 0.000001, 500, 1.0, 2.0, 0.5, 0.5, true);
						 
		alphaGlob = globalFit[0];
		betaGlob = globalFit[1];
	}
	
	const int rad = s;
	
	BufferedImage<float> astigFits(s,s,fc), astigNullFit(s,s,fc), astigData(s,s,fc);
	
	for (int f = 0; f < fc; f++)
	{
		double alpha, beta;
		
		if (global_astigmatism)
		{
			alpha = alphaGlob;
			beta = betaGlob;
		}
		else
		{
			std::vector<BufferedImage<float>> oneSpec(1);
			oneSpec[0] = avgTiles[f];
			
			AstigmatismFit afLocal(oneSpec, num_threads, r0_pix_ast, r1_pix_ast, Cs_px);
			
			std::vector<double> localFit = NelderMead::optimize(
				initial, afLocal, 0.01, 0.000001, 500, 1.0, 2.0, 0.5, 0.5, true);
							 
			alpha = localFit[0];
			beta = localFit[1];
		}
			
		if (diag)
		{
			std::vector<double> accum(rad), weight(rad);
			
			afGlob.averageAlongIso(f, alpha, beta, accum, weight, 2);
			
			std::stringstream sts;
			sts << f;
			
			BufferedImage<double> fit = afGlob.computeExpansion(f, alpha, beta, accum, 0);	
			
			NewStackHelper::insertSliceZ(
				Centering::fftwHalfToHumanFull(fit), astigFits, f);
			
			NewStackHelper::insertSliceZ(
				Centering::fftwHalfToHumanFull(avgTiles[f]), astigData, f);
			
			
			afGlob.averageAlongIso(f, 0.0, 0.0, accum, weight, 2);
			
			BufferedImage<double> nullFit = afGlob.computeExpansion(f, 0.0, 0.0, accum, 0);
			
			NewStackHelper::insertSliceZ(
				Centering::fftwHalfToHumanFull(nullFit), astigNullFit, f);
		}
	}
	
	std::string tag = global_astigmatism? "_global" : "_perFrame";
	
	astigData.write("debug/astigFitNM"+tag+"_data.mrc");
	astigFits.write("debug/astigFitNM"+tag+"_fits.mrc");
	astigNullFit.write("debug/astigFitNM"+tag+"_initial.mrc");
}

double TomoCtfFind::evaluate1D(double deltaF, int frame, double hand)
{
	const int s = tileSize;
	
	double cost = 0.0;
	
	for (int ty = 0; ty < th; ty++)
	for (int tx = 0; tx < tw; tx++)
	{		
		const double z = deltaF + hand * tileOffZ(tx, ty, frame);
				
		CTF ctf;
		ctf.setValues(z, z, 0.0, voltage, Cs, Q0, 0.0);
		
		for (int r = r0_pix; r <= r1_pix; r++)
		{
			const double ra = r / (double)(s * pixelSize);
			
			const double c = ctf.getCTF(ra, 0.0);
			const double d = radAvgNrm(tx, ty, frame)[r];
			
			cost -= frqWgh[r] * (c*c - 0.5) * d;
		}
	}
	
	return cost / (tw * th);
}

void TomoCtfFind::extractTiles()
{
	const int w = tiltSeries.xdim;
	const int h = tiltSeries.ydim;
	const int wh = w/2 + 1;
	
	const int s = tileSize;
	const int sh = s/2 + 1;
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<float> tempReal(s,s), tempRealBig(w,h);
		BufferedImage<tComplex<float>> tempComplex(sh,s), tempComplexBig(wh, h);
		
		tempRealBig = NewStackHelper::extractSliceZ(tiltSeries, f);
		FFT::FourierTransform(tempRealBig, tempComplexBig, FFT::Both);
		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < wh; x++)
		{
			const double xa = x / (double) (w * pixelSize);
			const double ya = (y < h/2? y : y - h) / (double) (h * pixelSize);
			
			const double rr = xa * xa + ya * ya;
			
			tempComplexBig(x,y) /= Damage::getWeight(frameDose[f], sqrt(rr)) + DAMAGE_EPS;
		}
		
		FFT::inverseFourierTransform(tempComplexBig, tempRealBig, FFT::Both);
		
		for (int ty = 0; ty < th; ty++)
		for (int tx = 0; tx < tw; tx++)
		{
			const int x0 = tx * tdx;
			const int y0 = ty * tdy;
			
			const int zz = f * tw * th + ty * tw + tx;
			
			for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				tempReal(x,y) = tempRealBig(x0 + x, y0 + y);
			}
			
			FFT::FourierTransform(tempReal, tempComplex, FFT::Both);
					
			for (int y = 0; y < s; y++)
			for (int x = 0; x < sh; x++)
			{
				tiles(x,y,zz) = tempComplex(x,y).norm();
			}
			
			const double xc = x0 + 0.5 * tileSize - w / 2.0;
			const double yc = y0 + 0.5 * tileSize - h / 2.0;
			
			d4Matrix A4 = projections[f];
			
			d2Matrix A2(
				A4(0,0), A4(0,1), 
				A4(1,0), A4(1,1));
			A2.invert();
			
			d2Vector a3(A4(2,0), A4(2,1));
			d2Vector ct(xc, yc);
			
			tileOffZ(tx, ty, f) = pixelSize * a3.dot(A2 * ct);
		}
	}
}

void TomoCtfFind::averageTiles()
{
	const int s = tileSize;
	const int sh = s/2 + 1;
	
	avgTiles.resize(fc);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		avgTiles[f] = BufferedImage<float>(sh,s);
		avgTiles[f].fill(0.f);
		
		for (int ty = 0; ty < th; ty++)
		for (int tx = 0; tx < tw; tx++)
		{
			const int zz = f * tw * th + ty * tw + tx;
			
			for (int y = 0; y < s; y++)
			for (int x = 0; x < sh; x++)
			{
				avgTiles[f](x,y) += tiles(x,y,zz);
			}
		}
		
		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			avgTiles[f](x,y) = avgTiles[f](x,y) / (tw * th);
		}
	}	
	
	avgTile = BufferedImage<float>(sh,s);
	avgTile.fill(0.f);
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < sh; x++)
	{
		double a = 0.0;
		
		for (int f = 0; f < fc; f++)
		{
			a += avgTiles[f](x,y);
		}
		
		avgTile(x,y) = a / fc;
	}
	
	if (diag)
	{
		avgTile.write("debug/avgTile.mrc");
	}
	
	double avgVal(0.0), avgWgh(0.0);
	
	const int b0 = 2, b1 = 1;
	
	for (int y = b0; y < s-b1; y++)
	for (int x = b0; x < sh; x++)
	{
		const double r = sqrt(x*x + y*y);
		
		if (r >= r0_pix && r <= r1_pix)
		{
			avgVal += avgTile(x,y);
			avgWgh += 1.0;
		}
	}
	
	avgVal /= avgWgh;
	
	avgTile /= avgVal;
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		avgTiles[f] /= avgVal;
		
		tiles.getSlabRef(f * tw * th, tw * th) /= avgVal;
	}
}

void TomoCtfFind::subtractBackground()
{
	const bool global_bg = true;
	const bool common_bg = true;
	
	const int s = tileSize;
	const int sh = s/2 + 1;
	
	const double sigmaK_Ang = 50.0;
	const double sigmaK_px = (tileSize * pixelSize) / sigmaK_Ang;
	const double tauK_px = sigmaK_px / 2.0;
	
	if (global_bg)
	{
		avgBackground = estimateBackground(avgTile, tauK_px, sigmaK_px, sh/3.0, diag);
		
		avgTileBgSub = avgTile - avgBackground;
		
		if (diag)
		{
			avgBackground.write("debug/global_bg.mrc");
			avgTileBgSub.write("debug/avgTile_bgsub.mrc");
		}
		
		#pragma omp parallel for num_threads(num_threads)	
		for (int f = 0; f < fc; f++)
		{
			// TODO: try figuring out the optimal scale per tile or per frame
			
			for (int ty = 0; ty < th; ty++)
			for (int tx = 0; tx < tw; tx++)
			{
				const int zz = f * tw * th + ty * tw + tx;
				RawImage<float> sl = tiles.getSliceRef(zz);
				
				sl -= avgBackground;
			}
		}
	}
	else
	{
		backgrounds.resize(fc);
		
		#pragma omp parallel for num_threads(num_threads)	
		for (int f = 0; f < fc; f++)
		{
			backgrounds[f] = estimateBackground(avgTiles[f], tauK_px, sigmaK_px, sh/3.0, f == 0);
			
			avgTiles[f] = avgTiles[f] - backgrounds[f];
			
			if (f == 0)
			{
				avgTiles[f].write("debug/avgTiles_f0.mrc");
				tiles.getSlabRef(f * tw * th, tw * th).write("debug/tiles_f0_0.mrc");
			}
			
					
			for (int ty = 0; ty < th; ty++)
			for (int tx = 0; tx < tw; tx++)
			{
				const int zz = f * tw * th + ty * tw + tx;
				RawImage<float> sl = tiles.getSliceRef(zz);
				
				if (common_bg)
				{
					sl -= backgrounds[f];
				}
				else
				{
					sl -= estimateBackground(sl, tauK_px, sigmaK_px, sh/3.0, false);
				}
			}
			
			if (f == 0)
			{
				tiles.getSlabRef(f * tw * th, tw * th).write("debug/tiles_f0_1.mrc");
			}
		}
	}
}

void TomoCtfFind::averageRadially()
{
	const int s = tileSize;
	const int sh = s/2 + 1;
				
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		for (int ty = 0; ty < th; ty++)
		for (int tx = 0; tx < tw; tx++)
		{
			const int zz = f * tw * th + ty * tw + tx;
			
			radAvg(tx, ty, f) = RadialAvg::fftwHalf_2D_lin(
						tiles.getSliceRef(zz), 2);
			
			radVar(tx, ty, f) = RadialAvg::variance_fftwHalf_2D_lin(
						tiles.getSliceRef(zz), radAvg(tx, ty, f), 2);
		}
	}
	
	std::vector<float> radAvgAll = RadialAvg::fftwHalf_2D_lin(avgTile, 2);	
	
	std::ofstream raDebug("debug/radAvgAll.dat");
	
	for (int r = 0; r < sh; r++)
	{
		raDebug << r << " " << radAvgAll[r] << "\n";
	}
		
	for (int f = 0; f < fc; f++)
	{
		double mu = 0.0, cwgh = 0.0;
		
		for (int ty = 0; ty < th; ty++)
		for (int tx = 0; tx < tw; tx++)	
		for (int r = r0_pix; r <= r1_pix; r++)
		{
			mu += frqWgh[r] * radAvg(tx, ty, f)[r];
			cwgh += frqWgh[r];
		}
		
		mu /= cwgh;
				
		double var = 0.0;
		
		for (int ty = 0; ty < th; ty++)
		for (int tx = 0; tx < tw; tx++)	
		for (int r = r0_pix; r <= r1_pix; r++)
		{
			double d = radAvg(tx, ty, f)[r] - mu;
			var += frqWgh[r] * d*d;
		}
		
		const double sd = sqrt(var / cwgh);
		
		for (int ty = 0; ty < th; ty++)
		for (int tx = 0; tx < tw; tx++)
		{
			radAvgNrm(tx, ty, f) = std::vector<float>(sh);
			
			for (int r = 0; r < sh; r++)
			{
				radAvgNrm(tx, ty, f)[r] = (radAvg(tx, ty, f)[r] - mu) / sd;
			}
		}
	}
	
	double mu_all(0.0), cwgh_all(0.0);
	
	for (int r = r0_pix; r <= r1_pix; r++)
	{
		mu_all += frqWgh[r] * radAvgAll[r];
		cwgh_all += frqWgh[r];
	}
	
	mu_all /= cwgh_all;
	
	double var_all = 0.0;
	
	for (int r = r0_pix; r <= r1_pix; r++)
	{
		double d = radAvgAll[r] - mu_all;
		var_all += frqWgh[r] * d*d;
	}
	
	const double sd_all = sqrt(var_all / cwgh_all);
	
	radAvgAllNrm.resize(sh);
	
	for (int r = 0; r < sh; r++)
	{
		radAvgAllNrm[r] = (radAvgAll[r] - mu_all) / sd_all;
	}
}

BufferedImage<float> TomoCtfFind::estimateBackground(
		const RawImage<float>& img_fftwHalf, 
		double sigma_blur_px, 
		double sigma_cent_px, 
		double sigma_out_px, 
		bool debug)
{
	const int s = img_fftwHalf.ydim;
	const int sh = img_fftwHalf.xdim;
	
	BufferedImage<float> out(sh,s);
	BufferedImage<float> weight(sh,s);
	
	const double sc2 = 2.0 * sigma_cent_px * sigma_cent_px;
	const double so2 = 2.0 * sigma_out_px * sigma_out_px;
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < sh; x++)
	{
		double xx = x;
		double yy = y < s/2? y : y - s;
		double rr = xx*xx + yy*yy;
		
		if (x < 2 || y < 2 || y > s-2)
		{
			out(x,y) = 0.0;
			weight(x,y) = 0.0;
		}
		else
		{
			const double env0 = 1.0 - exp(-rr / sc2);
			const double env1 = exp(-rr / so2);
					
			out(x,y) = env0 * env1 * img_fftwHalf(x,y);
			weight(x,y) = env0 * env1;
		}
	}
			
	BufferedImage<float> bgFull = Centering::fftwHalfToHumanFull(out);
	BufferedImage<float> wgFull = Centering::fftwHalfToHumanFull(weight);
	
	if (debug)
	{
		bgFull.write("debug/bgFull_f0_0.mrc");
		wgFull.write("debug/wgFull_f0_0.mrc");
	}
	
	bgFull = ImageFilter::Gauss2D(bgFull, 0, sigma_blur_px, true);
	wgFull = ImageFilter::Gauss2D(wgFull, 0, sigma_blur_px, true);
	
	if (debug)
	{
		bgFull.write("debug/bgFull_f0_1.mrc");
		wgFull.write("debug/wgFull_f0_1.mrc");
	}
	
	bgFull /= wgFull;
	
	if (debug)
	{
		bgFull.write("debug/bgFull_f0_2.mrc");
	}
	
	return Centering::humanFullToFftwHalf(bgFull);
}


