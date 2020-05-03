#include "tile_ctf_cost.h"
#include <src/jaz/image/centering.h>
#include <omp.h>

TileCtfCost::TileCtfCost(
	RawImage<float>& tiles, int tw, int th,
	RawImage<double>& tileOffZ,
	RawImage<float>& background, 
	double pixelSize, double voltage, double Cs, double Q0, 
	double hand,
	int num_threads, int r0, int r1)
	:	
		tw(tw), 
		th(th),
		hand(hand),
		num_threads(num_threads),
		tileOffZ(tileOffZ)
{
	tileCosts.resize(tw * th);
	
	const int bb = background.zdim == 1? 0 : 1;
	
	if (tiles.zdim < tw * th)
	{
		REPORT_ERROR_STR(
			"TileCtfCost::c'tor: insufficient number of tiles ("
			<< tiles.zdim << ".");
	}
	
	for (int ty = 0; ty < th; ty++)
	for (int tx = 0; tx < tw; tx++)
	{
		const int ti = ty * tw + tx;
		
		tileCosts[ti] = SpectralCtfCost(
			tiles.getSliceRef(ti), background.getSliceRef(bb * ti), 
			pixelSize, voltage, Cs, Q0, num_threads, r0, r1);
	}
}

double TileCtfCost::f(const std::vector<double> &x, void *tempStorage) const
{
	/*   layout of x:	 
	     [defocus, alpha, beta, rho, B] [zeta]_1 [zeta]_2 ... [zeta]_tc
	     [      0,     1,    2,   3, 4] [ 5 ]    [ 6 ]    ... [tc+4]
	*/
	
	const int tc = tw * th;
	
	std::vector<double> x_tile(6);
	
	const int pad = 1024;
	std::vector<double> out_par(pad * num_threads, 0.0);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int t = 0; t < tc; t++)
	{ 
		const SpectralCtfCost& scc = tileCosts[t];
		const int ti = omp_get_thread_num();
		
		// defocus is offset for each tile:
		x_tile[0] = scc.offsetDefocusParam(x[0], hand * tileOffZ[t]);
				
		// astigmatism and envelope are the same:
		for (int i = 1; i < 5; i++)
		{
			x_tile[i] = x[i];
		}
		
		// noise level is free for each tile: 
		x_tile[5] = x[5 + t];
		
		out_par[ti * pad] += scc.f(x_tile, tempStorage);
	}
	
	double out = 0.0;
	
	for (int t = 0; t < num_threads; t++)
	{
		out += out_par[t*pad];
	}
	
	return out;
}

void TileCtfCost::grad(
	const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	const int tc = tw * th;
	const int pc = tc + 5;
	
	const int tpc = 6;
	const int pad = 1024;
	
	std::vector<std::vector<double>> 
		x_tile(num_threads), grad_par_tile(num_threads), grad_par_out(num_threads);
	
	for (int t = 0; t < num_threads; t++)
	{ 
		x_tile[t] = std::vector<double>(tpc + pad, 0.0);
		grad_par_tile[t] = std::vector<double>(tpc + pad, 0.0);
		grad_par_out[t] = std::vector<double>(pc + pad, 0.0);
	}
	
	#pragma omp parallel for num_threads(num_threads)
	for (int t = 0; t < tc; t++)
	{ 
		const SpectralCtfCost& scc = tileCosts[t];
		
		const int ti = omp_get_thread_num();
		
		x_tile[ti][0] = scc.offsetDefocusParam(x[0], hand * tileOffZ[t]);
				
		for (int i = 1; i < 5; i++)
		{
			x_tile[ti][i] = x[i];
		}
		 
		x_tile[ti][5] = x[5 + t];
		
		scc.grad(x_tile[ti], grad_par_tile[ti], tempStorage);
		
		for (int i = 0; i < 5; i++)
		{
			grad_par_out[ti][i] += grad_par_tile[ti][i];
		}
		 
		grad_par_out[ti][5 + t] += grad_par_tile[ti][5];
	}
	
	for (int p = 0; p < pc; p++)
	{
		gradDest[p] = 0.0;
	}
	
	for (int ti = 0; ti < num_threads; ti++)
	{
		for (int p = 0; p < pc; p++)
		{
			gradDest[p] += grad_par_out[ti][p];
		}
	}
}

std::vector<double> TileCtfCost::getInitialParams(double defocus)
{
	const int tc = tw * th;
	
	SpectralCtfCost& scc0 = tileCosts[0];
	
	std::vector<double> in = scc0.getInitialParams(defocus);
	std::vector<double> out(5 + tc, in[6]);
	
	for (int i = 0; i < 5; i++)
	{
		out[i] = in[i];
	}
	
	return out;
}

std::vector<double> TileCtfCost::getInitialParams(const std::vector<double>& global)
{
	const int tc = tw * th;
	
	std::vector<double> out(5 + tc, global[6]);
	
	for (int i = 0; i < 5; i++)
	{
		out[i] = global[i];
	}
	
	return out;
}

BufferedImage<float> TileCtfCost::renderTable(const std::vector<double> &x)
{
	SpectralCtfCost& scc0 = tileCosts[0];
	
	const int tc = tw * th;
	const int s = scc0.spectrum.ydim;
	
	BufferedImage<float> out(s, s, tc);
	
	std::vector<double> x_tile(6);
	
	for (int t = 0; t < tc; t++)
	{
		SpectralCtfCost& scc = tileCosts[t];
		
		x_tile[0] = scc.offsetDefocusParam(x[0], hand * tileOffZ[t]);
		
		for (int i = 1; i < 5; i++)
		{
			x_tile[i] = x[i];
		}
		
		x_tile[5] = x[5 + t];
		
		BufferedImage<float> img0_t = scc.render(x_tile);
		BufferedImage<float> img_t = Centering::fftwHalfToHumanFull(img0_t);
				
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			out(x,y,t) = img_t(x,y);
		}
	}

	return out;
}

BufferedImage<float> TileCtfCost::renderCentralFit(const std::vector<double> &x)
{
	SpectralCtfCost& scc = tileCosts[0];
	
	BufferedImage<float> img0_t = scc.render(x);
	return Centering::fftwHalfToHumanFull(img0_t);
}

BufferedImage<float> TileCtfCost::renderAlignedSpectrum(const std::vector<double> &x)
{
	SpectralCtfCost& scc0 = tileCosts[0];
	
	const int tc = tw * th;
	const int s = scc0.spectrum.ydim;
	const int sh = scc0.spectrum.xdim;
	
	BufferedImage<double> accum(sh,s), weight(sh,s);
	
	std::vector<double> x_tile(6);
	
	for (int t = 0; t < tc; t++)
	{ 
		SpectralCtfCost& scc = tileCosts[t];
		
		x_tile[0] = scc.offsetDefocusParam(x[0], hand * tileOffZ[t]);
		
		for (int i = 1; i < 5; i++)
		{
			x_tile[i] = x[i];
		}
		
		x_tile[5] = x[5 + t];
		
		scc.addAlignedSpectrum(x[0], x_tile, accum, weight);
	}
	
	BufferedImage<float> div(sh,s);
			
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const double wg = weight(x,y);
		
		if (wg > 0.0)
		{
			div(x,y) = (float) (accum(x,y) / wg);
		}
		else 
		{
			div(x,y) = 0.f;
		}
	}
	
	return Centering::fftwHalfToHumanFull(div);
}

BufferedImage<float> TileCtfCost::renderAvgSpectrum()
{
	SpectralCtfCost& scc0 = tileCosts[0];
	
	const int tc = tw * th;
	const int s = scc0.spectrum.ydim;
	const int sh = scc0.spectrum.xdim;
	
	BufferedImage<float> avg(sh,s);
	avg.fill(0);
	
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		double a = 0.0;
		
		for (int t = 0; t < tc; t++)
		{ 
			SpectralCtfCost& scc = tileCosts[t];
			a += scc.spectrum(x,y);
		}
		
		avg(x,y) = (float) (a / tc);
	}
	
	return Centering::fftwHalfToHumanFull(avg);
}

std::vector<double> TileCtfCost::radialAverageAligned(const std::vector<double> &x)
{
	SpectralCtfCost& scc0 = tileCosts[0];
	
	const int tc = tw * th;
	const int s = scc0.spectrum.ydim;
	const int sh = scc0.spectrum.xdim;
	
	std::vector<double> accum(sh,s), weight(sh,s);
	
	std::vector<double> x_tile(6);
	
	for (int t = 0; t < tc; t++)
	{ 
		SpectralCtfCost& scc = tileCosts[t];
		
		x_tile[0] = scc.offsetDefocusParam(x[0], hand * tileOffZ[t]);
		
		for (int i = 1; i < 5; i++)
		{
			x_tile[i] = x[i];
		}
		
		x_tile[5] = x[5 + t];
		
		scc.addAlignedRadialAverage(x[0], x_tile, accum, weight);
	}
	
	std::vector<double> div(sh,0.0);
			
	for (int r = 0; r < sh; r++)
	{
		const double wg = weight[r];
		
		if (wg > 0.0)
		{
			div[r] = accum[r] / wg;
		}
	}
	
	return div;	
}

std::vector<double> TileCtfCost::radialAverage(const std::vector<double> &x)
{
	SpectralCtfCost& scc0 = tileCosts[0];
	
	const int tc = tw * th;
	const int s = scc0.spectrum.ydim;
	const int sh = scc0.spectrum.xdim;
	
	std::vector<double> accum(sh,s), weight(sh,s);
	
	std::vector<double> x_tile(6);
	
	for (int t = 0; t < tc; t++)
	{ 
		SpectralCtfCost& scc = tileCosts[t];
		
		x_tile[0] = scc.offsetDefocusParam(x[0], hand * tileOffZ[t]);
		
		for (int i = 1; i < 5; i++)
		{
			x_tile[i] = x[i];
		}
		
		x_tile[5] = x[5 + t];
		
		scc.addRadialAverage(accum, weight);
	}
	
	std::vector<double> div(sh,0.0);
			
	for (int r = 0; r < sh; r++)
	{
		const double wg = weight[r];
		
		if (wg > 0.0)
		{
			div[r] = accum[r] / wg;
		}
	}
	
	return div;		
}


std::vector<double> TileCtfCost::radialAverageFit(const std::vector<double> &x)
{
	SpectralCtfCost& scc0 = tileCosts[0];
	
	const int tc = tw * th;
	const int s = scc0.spectrum.ydim;
	const int sh = scc0.spectrum.xdim;
	
	std::vector<double> accum(sh,s), weight(sh,s);
	
	std::vector<double> x_tile(6);
	
	for (int t = 0; t < tc; t++)
	{ 
		SpectralCtfCost& scc = tileCosts[t];
		
		x_tile[0] = x[0]; // not ofsetting Z
		
		for (int i = 1; i < 5; i++)
		{
			x_tile[i] = x[i];
		}
		
		x_tile[5] = x[5 + t];
		
		scc.addRadialAverageFit(x_tile, accum, weight);
	}
	
	std::vector<double> div(sh,0.0);
			
	for (int r = 0; r < sh; r++)
	{
		const double wg = weight[r];
		
		if (wg > 0.0)
		{
			div[r] = accum[r] / wg;
		}
	}
	
	return div;		
}
