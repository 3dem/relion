#include "tilt_space_blob_fit.h"
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/normalization.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>

using namespace gravis;
using namespace au::edu::anu::qm::ro;


TiltSpaceBlobFit::TiltSpaceBlobFit(
	int sh_bands, 
	double lambda,
	const RawImage<float>& correlation,
	const RawImage<d3Vector>& directions_xz)
	:
	sh_bands(sh_bands),
	lambda(lambda),
	correlation(correlation)
{
	SphericalHarmonics* SH = new SphericalHarmonics(sh_bands);
	
	const int w  = correlation.xdim;
	const int fc = correlation.zdim;
	const int SH_params = (sh_bands + 1) * (sh_bands + 1);
	
	basis.resize(SH_params, w, fc);
	
	for (int f = 0; f < fc; f++)
	for (int x = 0; x < w; x++)
	{
		const d3Vector v = directions_xz(x,0,f);
		const double phi = atan2(v.y, v.x);
		
		SH->computeY(sh_bands, v.z, phi, &basis(0,x,f));
	}
}

double TiltSpaceBlobFit::f(
        const std::vector<double>& x, 
        void* tempStorage) const
{
	const int w  = correlation.xdim;
	const int h  = correlation.ydim;
	const int fc = correlation.zdim;
	const int bc = basis.xdim;
	
	double corr_sum = 0.0;
	
	for (int f = 0; f < fc; f++)
	for (int xx = 0; xx < w; xx++)
	{
		double Y = 0.0;
		
		for (int b = 0; b < bc; b++)
		{
			Y += x[b] * basis(b,xx,f);
		}
		
		int h0 = (int) std::floor(Y);
		int h1 = h0 + 1;
		
		const double dh = Y - h0;
		
		if (h0 < 0) h0 = 0;
		else if (h0 > h - 1) h0 = h - 1;
		
		if (h1 < 0) h1 = 0;
		else if (h1 > h - 1) h1 = h - 1;
		
		const float cc0 = correlation(xx, h0, f);
		const float cc1 = correlation(xx, h1, f);
		
		corr_sum += (1 - dh) * cc0 + dh * cc1;
	}
	
	double out = -corr_sum / (fc * (double)w);
	
	for (int b = 1; b < bc; b++)
	{
		out += lambda * std::floor(sqrt((double)b)) * x[b] * x[b];
	}
	
	return out;
}

void TiltSpaceBlobFit::grad(
		const std::vector<double>& x,
		std::vector<double>& gradDest,
		void *tempStorage) const
{
	const int w  = correlation.xdim;
	const int h  = correlation.ydim;
	const int fc = correlation.zdim;
	const int bc = basis.xdim;
	
	for (int b = 0; b < bc; b++)
	{
		gradDest[b] = 0.0;
	}
	
	for (int f = 0; f < fc; f++)
	for (int xx = 0; xx < w; xx++)
	{
		double Y = 0.0;
		
		for (int b = 0; b < bc; b++)
		{
			Y += x[b] * basis(b,xx,f);
		}
		
		int h0 = (int) std::floor(Y);
		int h1 = h0 + 1;
		
		if (h0 < 0) h0 = 0;
		else if (h0 > h - 1) h0 = h - 1;
		
		if (h1 < 0) h1 = 0;
		else if (h1 > h - 1) h1 = h - 1;
		
		const float cc0 = correlation(xx, h0, f);
		const float cc1 = correlation(xx, h1, f);
		const float dcc_dY = cc1 - cc0;
		
		for (int b = 0; b < bc; b++)
		{
			gradDest[b] -= dcc_dY * basis(b,xx,f)  / (fc * (double)w);
		}
	}
	
	for (int b = 1; b < bc; b++)
	{
		gradDest[b] += 2.0 * lambda * std::floor(sqrt((double)b)) * x[b];
	}
}

double TiltSpaceBlobFit::estimateInitialHeight()
{
	const int w  = correlation.xdim;
	const int h  = correlation.ydim;
	const int fc = correlation.zdim;
	
	std::vector<double> sum_per_y(h, 0.0);
	
	for (int f = 0;  f < fc; f++)
	for (int y = 0;  y < h;  y++)
	for (int xx = 0; xx < w; xx++)
	{
		sum_per_y[y] += correlation(xx,y,f);
	}
	
	int best_y = h/2;
	double max_sum = 0.0;
	
	for (int y = 0; y < h; y++)
	{
		if (sum_per_y[y] > max_sum)
		{
			max_sum = sum_per_y[y];
			best_y = y;
		}
	}
	
	return best_y;
}

BufferedImage<float> TiltSpaceBlobFit::drawSolution(
        const std::vector<double> &x,
        const RawImage<float>& map)
{
	const int w  = correlation.xdim;
	const int h  = correlation.ydim;
	const int fc = correlation.zdim;
	const int bc = basis.xdim;
	
	BufferedImage<float> out = map;
	const float mu = Normalization::computeMean(map);
	const float var = Normalization::computeVariance(map, mu);
	const float val = mu + 6 * sqrt(var);
	
	for (int f = 0; f < fc; f++)
	for (int xx = 0; xx < w; xx++)
	{
		double Y = 0.0;
		
		for (int b = 0; b < bc; b++)
		{
			Y += x[b] * basis(b,xx,f);
		}
		
		int h0 = (int) std::floor(Y);
		int h1 = h0 + 1;
		
		const double dh = Y - h0;
		
		if (h0 >= 0 && h0 < h)
		{
			out(xx,h0,f) += (float)(1.0 - dh) * val;
		}
		
		if (h1 >= 0 && h1 < h)
		{
			out(xx,h1,f) += (float) dh * val;
		}
	}
	
	return out;
}

int TiltSpaceBlobFit::getParameterCount()
{
	return basis.xdim;
}

BufferedImage<float> TiltSpaceBlobFit::computeTiltSpaceMap(
        d3Vector sphere_position, 
        double mean_radius_full, 
		double radius_range,
        double binning, 
        const RawImage<float>& preweighted_stack, 
        const std::vector<d4Matrix>& projections)
{
	const int fc = projections.size();
		
	const double perimeter_length_full = 2 * PI * mean_radius_full;
	
	int w_map = perimeter_length_full / binning;
	int h_map = radius_range / binning;
	
	w_map += w_map % 2;
	h_map += h_map % 2;
	
	const int w_stack = preweighted_stack.xdim;
	const int h_stack = preweighted_stack.ydim;
	
	const double min_radius_full = mean_radius_full - radius_range/2;
	
	BufferedImage<float> map(w_map, h_map, fc);
	
	for (int z = 0; z < fc; z++)
	for (int y = 0; y < h_map; y++)
	for (int x = 0; x < w_map; x++)
	{
		const int f0 = z;
		const d4Matrix& A = projections[f0];
		
		d3Vector dir_x(A(0,0), A(0,1), A(0,2));
		d3Vector dir_y(A(1,0), A(1,1), A(1,2));
		
		dir_x.normalize();
		dir_y.normalize();
		
		const double phi = 2 * PI * x / (double) w_map;
		const double r = min_radius_full + radius_range * y / (double) h_map;
		
		const d3Vector pos = sphere_position + r * (cos(phi) * dir_x + sin(phi) * dir_y);
		
		float sum = 0.f;
		float weight = 0.f;
		
		for (int f = 0; f < fc; f++)
		{
			const d4Vector pi = projections[f] * d4Vector(pos);

			if (pi.x >= 0.0 && pi.x < w_stack && pi.y >= 0.0 && pi.y < h_stack)
			{
				sum += Interpolation::linearXY_clip(preweighted_stack, pi.x, pi.y, f);
				weight += 1.0;
			}
		}
		
		if (weight > 0)
		{
			map(x,y,z) = sum / weight;
		}
		else
		{
			map(x,y,z) = 0;
		}     
	}
	
	return map;
}

BufferedImage<d3Vector> TiltSpaceBlobFit::computeDirectionsXZ(
        double mean_radius_full, 
        double binning, 
        const std::vector<d4Matrix>& projections)
{
	const int fc = projections.size();
		
	const double perimeter_length_full = 2 * PI * mean_radius_full;
	
	int w_map = perimeter_length_full / binning;
	
	w_map += w_map % 2;
	
	BufferedImage<d3Vector> directions_XZ(w_map, 1, fc);
	
	for (int f = 0; f < fc; f++)
	for (int x = 0; x < w_map; x++)
	{
		const int f0 = f;
		const d4Matrix& A = projections[f0];
		
		d3Vector dir_x(A(0,0), A(0,1), A(0,2));
		d3Vector dir_y(A(1,0), A(1,1), A(1,2));
		
		dir_x.normalize();
		dir_y.normalize();
		
		const double phi = 2 * PI * x / (double) w_map;
		
		directions_XZ(x,0,f) = cos(phi) * dir_x + sin(phi) * dir_y;  
	}
	
	return directions_XZ;
}

BufferedImage<float> TiltSpaceBlobFit::visualiseBlob(
        const std::vector<double> parameters, 
        double mean_radius_full, 
        double radius_range, 
        double binning, 
        const RawImage<float>& preweighted_stack, 
        const std::vector<d4Matrix>& projections)
{
	const int fc = projections.size();
		
	const double perimeter_length_full = 2 * PI * mean_radius_full;
	
	const d3Vector centre(parameters[0], parameters[1], parameters[2]);
	
	const int SH_params = parameters.size() - 3;
	
	const int SH_bands = (int) sqrt(SH_params - 1);
	
	const float mu = Normalization::computeMean(preweighted_stack);
	const float var = Normalization::computeVariance(preweighted_stack, mu);
	const float draw_val = mu + 6 * sqrt(var);
	
	SphericalHarmonics SH(SH_bands);
	
	int w_map = perimeter_length_full / binning;
	int h_map = radius_range / binning;
	
	w_map += w_map % 2;
	h_map += h_map % 2;
	
	const int w_stack = preweighted_stack.xdim;
	const int h_stack = preweighted_stack.ydim;
		
	BufferedImage<float> map(w_map, h_map, fc);
	std::vector<double> Y(SH_params);
	
	const int y0 = h_map / 2;
	
	for (int z = 0; z < fc; z++)
	for (int y = 0; y < h_map; y++)
	for (int x = 0; x < w_map; x++)
	{
		const int f0 = z;
		const d4Matrix& A = projections[f0];
		
		d3Vector dir_x(A(0,0), A(0,1), A(0,2));
		d3Vector dir_y(A(1,0), A(1,1), A(1,2));
		
		dir_x.normalize();
		dir_y.normalize();
		
		const double phi_2D = 2 * PI * x / (double) w_map;
		const d3Vector d = cos(phi_2D) * dir_x + sin(phi_2D) * dir_y;
		        
		SH.computeY(SH_bands, d.z, atan2(d.y, d.x), &Y[0]);
		
		double rad = 0.0;
		
		for (int i = 0; i < SH_params; i++)
		{
			rad += parameters[i+3] * Y[i];
		}
		        
		const d3Vector pos = centre + (rad + radius_range * ((y - y0)/ (double) h_map)) * d;
		
		        
		float sum = 0.f;
		float weight = 0.f;
		
		for (int f = 0; f < fc; f++)
		{
			const d4Vector pi = projections[f] * d4Vector(pos);

			if (pi.x >= 0.0 && pi.x < w_stack && pi.y >= 0.0 && pi.y < h_stack)
			{
				sum += Interpolation::linearXY_clip(preweighted_stack, pi.x, pi.y, f);
				weight += 1.0;
			}
		}
		
		if (weight > 0)
		{
			map(x,y,z) = sum / weight;
		}
		else
		{
			map(x,y,z) = 0;
		}
		
		if (y == y0)
		{
			map(x,y,z) = draw_val;
		}
	}
	
	return map;
}

