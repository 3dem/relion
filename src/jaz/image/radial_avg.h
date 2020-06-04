#ifndef RADIAL_AVG_H
#define RADIAL_AVG_H

#include <vector>
#include <src/image.h>
#include "buffered_image.h"

class RadialAvg
{
	public:

		template<typename T>
		static std::vector<T> fftwHalf_2D_lin(const RawImage<T>& img, int boundary = 0);

		template<typename T>
		static std::vector<T> fftwHalf_3D_lin(const RawImage<T>& img);

		template<typename T>
		static T interpolate_FftwHalf_3D_lin(
				double x, double y, double z,
				int w, int h, int d,
				const std::vector<T>& avg);

		static double get1DIndex(
				double x, double y, double z,
				int w, int h, int d);



		template<typename T>
		static std::vector<T> variance_fftwHalf_2D_lin(
				const RawImage<T>& img,
				const std::vector<T>& avg,
				int boundary = 0);

		template<typename T>
		static std::vector<T> fftwHalf_2D_NN(const RawImage<T>& img, int boundary = 0);

		template<typename T>
		static void toFftwHalf_2D_lin(const std::vector<T>& vec, int wh, int h, RawImage<T>& dest, double scale = 1.0);

		template<typename T>
		static BufferedImage<T> toFftwHalf_2D_lin(const std::vector<T>& vec, int wh, int h, double scale = 1.0);

		template<typename T>
		static std::vector<std::pair<T,T>> radialAverageAndStdDevFFTW_2D(const Image<T>& map);


		template<typename T>
		static std::vector<T> radialAverageFFTW_2D(const Image<T>& map);

		template<typename T>
		static std::vector<T> radialSumFFTW_2D(const Image<T>& map);
};

template<typename T>
std::vector<T> RadialAvg::fftwHalf_2D_lin(const RawImage<T>& img, int boundary)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int b = w;

	std::vector<T> avg(b, 0.0);
	std::vector<double> wgh(b, 0.0);

	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;

	for (int y = b0; y < h-b1; y++)
	for (int x = b0; x < w-b1; x++)
	{
		const double xx = x;
		const double yy = (y >= h/2)? y - h: y;

		const double r = sqrt(xx * xx + yy * yy);

		const int r0 = (int) r;
		const int r1 = r0 + 1;

		const T val = img(x,y);

		const double w1 = r - r0;
		const double w0 = 1.0 - w1;

		if (r0 < b)
		{
			avg[r0] += w0 * val;
			wgh[r0] += w0;
		}

		if (r1 < b)
		{
			avg[r1] += w1 * val;
			wgh[r1] += w1;
		}
	}

	for (int i = 0; i < b; i++)
	{
		if (wgh[i] > 0.0) avg[i] /= wgh[i];
	}

	return avg;
}

template<typename T>
std::vector<T> RadialAvg::fftwHalf_3D_lin(const RawImage<T>& img)
{
	const int wh = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;

	std::vector<T> avg(wh, 0.0);
	std::vector<double> wgh(wh, 0.0);

	for (int z = 0; z < d;  z++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x;
		const double yy = (y >= h/2)? y - h: y;
		const double zz = (z >= d/2)? z - d: z;

		const double r = sqrt(xx * xx + yy * yy + zz * zz);

		const int r0 = (int) r;
		const int r1 = r0 + 1;

		const T val = img(x,y);

		const double w1 = r - r0;
		const double w0 = 1.0 - w1;

		if (r0 < wh)
		{
			avg[r0] += w0 * val;
			wgh[r0] += w0;
		}

		if (r1 < wh)
		{
			avg[r1] += w1 * val;
			wgh[r1] += w1;
		}
	}

	for (int i = 0; i < wh; i++)
	{
		if (wgh[i] > 0.0) avg[i] /= wgh[i];
	}

	return avg;
}

template<typename T>
T RadialAvg::interpolate_FftwHalf_3D_lin(
		double x, double y, double z,
		int w, int h, int d,
		const std::vector<T>& avg)
{
	const int s = avg.size();

	const double xx = x;
	const double yy = (y >= h/2)? y - h: y;
	const double zz = (z >= d/2)? z - d: z;

	const double r = sqrt(xx * xx + yy * yy + zz * zz);

	const int r0 = (int) r;
	const int r1 = r0 + 1;

	if (r1 >= s) return avg[s-1];

	const T v1 = avg[r1];
	const T v0 = avg[r0];

	const double f = r - r0;

	return f * v1 + (1.0 - f) * v0;
}

template<typename T>
std::vector<T> RadialAvg::variance_fftwHalf_2D_lin(
		const RawImage<T>& img,
		const std::vector<T>& avg,
		int boundary)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int b = w;

	std::vector<T> var(b, 0.0);
	std::vector<double> wgh(b, 0.0);

	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;

	for (int y = b0; y < h-b1; y++)
	for (int x = b0; x < w-b1; x++)
	{
		const double xx = x;
		const double yy = (y >= h/2)? y - h: y;

		const double r = sqrt(xx * xx + yy * yy);

		const int r0 = (int) r;
		const int r1 = r0 + 1;

		if (r1 >= b) continue;


		const double w1 = r - r0;
		const double w0 = 1.0 - w1;

		const T av = w0 * avg[r0] + w1 * avg[r1];

		const T diff = img(x,y) - av;
		const T d2 = diff * diff;

		var[r0] += w0 * d2;
		wgh[r0] += w0;

		var[r1] += w1 * d2;
		wgh[r1] += w1;
	}

	for (int i = 0; i < b; i++)
	{
		if (wgh[i] > 1.0) var[i] /= (wgh[i] - 1.0);
	}

	return var;
}

template<typename T>
std::vector<T> RadialAvg::fftwHalf_2D_NN(const RawImage<T>& img, int boundary)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int b = w;

	std::vector<T> avg(b, 0.0);
	std::vector<T> wgh(b, 0.0);

	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;

	for (int yy = b0; yy < h-b1; yy++)
	for (int xx = b0; xx < w-b1; xx++)
	{
		const double x = xx;
		const double y = yy < h/2.0? yy : yy - h;
		const double rd = sqrt(x*x + y*y);

		const int r = (int)(rd + 0.5);
		const double ct = (xx == 0)? 1.0 : 2.0;

		if (r < b)
		{
			avg[r] += ct * img(xx,yy);
			wgh[r] += ct;
		}
	}

	for (int i = 0; i < b; i++)
	{
		if (wgh[i] > 0.0)
		{
			avg[i] /= wgh[i];
		}
	}

	return avg;
}


template<typename T>
void RadialAvg::toFftwHalf_2D_lin(const std::vector<T>& vec, int wh, int h, RawImage<T>& dest, double scale)
{
	const int b = vec.size();

	for (int yy = 0; yy < h; yy++)
	for (int xx = 0; xx < wh; xx++)
	{
		const double x = xx;
		const double y = yy < h/2.0? yy : yy - h;
		const double rd = scale * sqrt(x*x + y*y);

		const int r0 = (int)(rd);
		const int r1 = r0 + 1;

		if (r1 >= b)
		{
			dest(xx,yy) = vec[b-1];
		}
		else
		{
			const double t = rd - r0;

			dest(xx,yy) = (1 - t) * vec[r0] + t * vec[r1];
		}
	}
}

template<typename T>
BufferedImage<T> RadialAvg::toFftwHalf_2D_lin(const std::vector<T>& vec, int wh, int h, double scale)
{
	BufferedImage<T> out(wh,h);
	out.fill(T(0));

	toFftwHalf_2D_lin(vec, wh, h, out, scale);

	return out;
}

template<typename T>
std::vector<std::pair<T,T>> RadialAvg::radialAverageAndStdDevFFTW_2D(const Image<T>& map)
{
	const int w = map.data.xdim;
	const int h = map.data.ydim;
	const int b = w;

	std::vector<T> avg(b, 0.0);
	std::vector<T> wgh(b, 0.0);
	std::vector<T> var(b, 0.0);

	for (int yy = 0; yy < h; yy++)
	for (int xx = 0; xx < w; xx++)
	{
		double x = xx;
		double y = yy < h/2.0? yy : yy - h;
		double rd = sqrt(x*x + y*y);

		int r = (int)(rd+0.5);

		if (r < b)
		{
			avg[r] += DIRECT_A2D_ELEM(map.data, yy, xx);
			wgh[r] += 1.0;
		}
	}

	for (int i = 0; i < b; i++)
	{
		if (wgh[i] > 0.0)
		{
			avg[i] /= wgh[i];
		}
	}

	for (int yy = 0; yy < h; yy++)
	for (int xx = 0; xx < w; xx++)
	{
		double x = xx;
		double y = yy < h/2.0? yy : yy - h;
		double rd = sqrt(x*x + y*y);

		int r = (int)(rd+0.5);

		T mu = avg[r];
		T v = DIRECT_A2D_ELEM(map.data, yy, xx) - mu;

		if (r < b)
		{
			var[r] += v*v;
		}
	}

	for (int i = 0; i < b; i++)
	{
		if (wgh[i] > 1.0)
		{
			var[i] /= (wgh[i]-1);
		}
	}

	std::vector<std::pair<T,T>> out(b);

	for (int i = 0; i < b; i++)
	{
		out[i] = std::make_pair(avg[i], sqrt(var[i]));
	}

	return out;
}

template<typename T>
std::vector<T> RadialAvg::radialAverageFFTW_2D(const Image<T>& map)
{
	const int w = map.data.xdim;
	const int h = map.data.ydim;
	const int b = w;

	std::vector<T> avg(b, 0.0);
	std::vector<T> wgh(b, 0.0);

	for (int yy = 0; yy < h; yy++)
	for (int xx = 0; xx < w; xx++)
	{
		double x = xx;
		double y = yy < h/2.0? yy : yy - h;
		double rd = sqrt(x*x + y*y);

		int r = (int)(rd+0.5);

		if (r < b)
		{
			avg[r] += DIRECT_A2D_ELEM(map.data, yy, xx);
			wgh[r] += 1.0;
		}
	}

	for (int i = 0; i < b; i++)
	{
		if (wgh[i] > 0.0)
		{
			avg[i] /= wgh[i];
		}
	}

	return avg;
}

template<typename T>
std::vector<T> RadialAvg::radialSumFFTW_2D(const Image<T>& map)
{
	const int w = map.data.xdim;
	const int h = map.data.ydim;
	const int b = w;

	std::vector<T> sum(b, 0.0);

	for (int yy = 0; yy < h; yy++)
	for (int xx = 0; xx < w; xx++)
	{
		double x = xx;
		double y = yy < h/2.0? yy : yy - h;
		double rd = sqrt(x*x + y*y);

		int r = (int)(rd+0.5);

		if (r < b)
		{
			sum[r] += DIRECT_A2D_ELEM(map.data, yy, xx);
		}
	}

	return sum;
}

#endif
