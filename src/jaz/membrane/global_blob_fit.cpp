#include "global_blob_fit.h"
#include "membrane_segmentation.h"
#include <src/jaz/math/fft.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/util/drawing.h>

using namespace gravis;


std::pair<double,std::vector<double>> GlobalBlobFit2D::fit(
        const std::vector<double> &initial_parameters, 
        double initial_radius, 
        int max_frequencies,
        const RawImage<float>& image, 
        const RawImage<float>& mask)
{
	const bool debug = false;
	
	std::vector<double> parameters(initial_parameters);
	
	for (int i = 2; i < parameters.size(); i++)
	{
		parameters[i] = 0.0;
	}
	
	Blob2D blob(parameters, initial_radius/2.0);
	
	std::pair<BufferedImage<float>, BufferedImage<float>> polar = 
	        blob.transformToPolar(image, mask, 2*initial_radius);
	
	BufferedImage<float> polar_image = polar.first;
	BufferedImage<float> polar_mask = polar.second;
	
	if (debug)
	{
		polar_image.write("DEBUG_glob_polar_image.mrc");
		polar_mask.write("DEBUG_glob_polar_mask.mrc");
	}
		
	const int w = polar_image.xdim;
	const int wh = w/2 + 1;
	const int h = polar_image.ydim;
	const int falloff = 100; // falloff
	const int width = 10; // width
	const double spacing = 33.0; // spacing
	const double ratio = 5.0; // ratio
	const double margin_out = h / 6.0;
	const double margin_in = h / 3.0;
		
	const double step_cost = 15;
	const int max_step = 5;
	const int num_freqs = max_frequencies;
		
	BufferedImage<float> kernel = MembraneSegmentation::constructMembraneKernel(
		w, h, 1, falloff, width, spacing, ratio, 0.0);   
	 
	if (debug)
	{
		kernel.write("DEBUG_glob_kernel.mrc");
	}
	
	
	BufferedImage<fComplex> kernel_FS, polar_image_FS, polar_mask_FS, correlation_FS;
	
	FFT::FourierTransform(kernel, kernel_FS, FFT::Both);
	FFT::FourierTransform(polar_image, polar_image_FS, FFT::Both);
	FFT::FourierTransform(polar_mask, polar_mask_FS, FFT::Both);
	
	correlation_FS.resize(wh,h);
	
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		correlation_FS(x,y) = polar_image_FS(x,y) * kernel_FS(x,y).conj();
	}
	
	BufferedImage<float> correlation;
	
	FFT::inverseFourierTransform(correlation_FS, correlation, FFT::Both);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double yy = y < h/2? y : y - h;
		
		if (yy > 0)
		{
			correlation(x,y) *= (1.0 - exp(-yy*yy / (2 * margin_in * margin_in)));
		}
		else
		{
			correlation(x,y) *= (1.0 - exp(-yy*yy / (2 * margin_out * margin_out)));
		}		
	}
	
	const float var = Normalization::computeVariance(correlation, 0.f);
	correlation /= sqrt(var);
	
	if (debug)
	{
		correlation.write("DEBUG_glob_correlation.mrc");
	}
	
	BufferedImage<double> min_cost(w,h);
	BufferedImage<int> best_origin(w,h);
	
	for (int y = 0; y < h; y++)
	{
		min_cost(0,y) = -correlation(0,y);
	}
	
	for (int x = 1; x < w; x++)
	{
		for (int y = 0; y < h; y++)
		{
			double min_hyp_cost = std::numeric_limits<double>::max();
			int best_hyp_origin = y;
			
			for (int dy = -max_step; dy <= max_step; dy++)
			{
				const int yy = y + dy;
				
				if (yy >= 0 && yy < h)
				{
					//const double hyp_cost = min_cost((x + w - 1) % w, yy) + std::abs(dy) * step_cost;
					const double hyp_cost = min_cost((x + w - 1) % w, yy) + dy * dy * step_cost;
					
					if (hyp_cost < min_hyp_cost) 
					{
						min_hyp_cost = hyp_cost;
						best_hyp_origin = yy;
					}
				}
				
				min_cost(x,y) = min_hyp_cost - correlation(x, y);
				best_origin(x,y) = best_hyp_origin;
			}
		}
	}
	
	// determine y(0) of all y(w-1)
	// add the cost of closing up from y(w-1) to y(0) to the cost of y(w-1)
		
	for (int y = 0; y < h; y++)
	{
		int yy = best_origin(w-1,y);
		
		for (int x = w-2; x > 0; x--)
		{
			yy = best_origin(x,yy);
		}
		
		const int y0 = yy;
		
		//min_cost(w-1, y) += std::abs(y0-y) * step_cost;
		min_cost(w-1, y) += (y0-y) * (y0-y) * step_cost;
	}
	
	
	double min_total_cost = std::numeric_limits<double>::max();
	int best_y = 0;
	        
	for (int y = 0; y < h; y++)
	{	
		if (min_cost(w-1, y) < min_total_cost)
		{
			min_total_cost = min_cost(w-1, y);
			best_y = y;
		}
	}
	
	if (debug)
	{
		min_cost.write("DEBUG_glob_min_cost.mrc");
	}
	
	
	std::vector<int> optimal_y(w);
	
	int current_y = best_y;
	optimal_y[w-1] = current_y;
	
	for (int x = w-2; x >= 0; x--)
	{
		current_y = best_origin(x+1, current_y);
		optimal_y[x] = current_y;		
	}
	
	
	double mean_radius = 0.0;
	
	for (int x = 0; x < w; x++)
	{
		mean_radius += optimal_y[x];
	}
	
	mean_radius /= (double) w;
	
	
	std::vector<fComplex> amplitudes(num_freqs+1);
	
	
	for (int f = 0; f < num_freqs+1; f++)
	{
		const double N = 1 + f;
		
		double dot_sin = 0.0;
		double dot_cos = 0.0;
		double sin_2 = 0.0;
		double cos_2 = 0.0;
		
		for (int x = 0; x < w; x++)
		{
			const double phi = 2 * PI * N * x / (double) w;
			
			dot_sin += sin(phi) * (optimal_y[x] - mean_radius);
			dot_cos += cos(phi) * (optimal_y[x] - mean_radius);
			
			sin_2 += sin(phi) * sin(phi);
			cos_2 += cos(phi) * cos(phi);
		}
		
		amplitudes[f].real = dot_cos / cos_2;
		amplitudes[f].imag = dot_sin / sin_2;
	}
	
	if (debug)
	{
		for (int x = 0; x < w; x++)
		{
			double radius = mean_radius;
			
			for (int f = 0; f < num_freqs+1; f++)
			{
				const double N = 1 + f;
				const double phi = 2 * PI * N * x / (double) w;
				
				radius += 
					amplitudes[f].real * cos(phi) +
					amplitudes[f].imag * sin(phi);
			}
			
			const int ri = (int) radius;
			const float rf = radius - ri;
			
			if (ri >= 0 && ri < h)
			{
				correlation(x, ri) += 15 * (1 - rf);
			}
			
			if (ri >= -1 && ri < h-1)
			{
				correlation(x, ri+1) += 15 * rf;
			}
			
			if ((x / 4 + 2) % 2)
			{
				correlation(x, optimal_y[x]) = 12;
			}
		}
		
		correlation.write("DEBUG_glob_solution.mrc");
	}
	
	
	std::vector<double> out(2 * num_freqs + 2);
	
	out[0] = initial_parameters[0] + amplitudes[0].real;
	out[1] = initial_parameters[1] + amplitudes[0].imag;
	
	for (int f = 0; f < num_freqs; f++)
	{
		out[2 * f + 2] = -amplitudes[f+1].real;
		out[2 * f + 3] = -amplitudes[f+1].imag;
	}
	
	if (debug)
	{	
		BufferedImage<float> test_image = drawOutline(out, mean_radius, image);		
		test_image.write("DEBUG_glob_outline.mrc");
	}
	     
	return std::make_pair(mean_radius, out);
}

BufferedImage<float> GlobalBlobFit2D::drawOutline(
        const std::vector<double> parameters,
        double mean_radius,
        const BufferedImage<float>& image)
{	
	const int phi_count = 2 * PI * mean_radius;
	
	DelineatedBlob2D test_blob(Blob2D(parameters, 0), mean_radius);	
	BufferedImage<float> test_image = image;
	
	for (int x = 0; x < phi_count; x++)
	{
		const double phi = 2 * PI * x / (double) phi_count;
		
		const d2Vector p = test_blob.getOutlinePoint(phi);
		
		Drawing::drawPoint(p, 10.f, test_image);	
	}
	
	return test_image;
}
	
