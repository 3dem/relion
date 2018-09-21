#ifndef RADIAL_AVG_H
#define RADIAL_AVG_H

#include <vector>
#include <src/image.h>

class RadialAvg
{
	public:
		
		template<typename T>
		static std::vector<std::pair<T,T>> radialAverageAndStdDevFFTW_2D(const Image<T>& map)
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
		static std::vector<T> radialAverageFFTW_2D(const Image<T>& map)
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
		static std::vector<T> radialSumFFTW_2D(const Image<T>& map)
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
};

#endif
