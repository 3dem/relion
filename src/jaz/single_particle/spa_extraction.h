#ifndef SPA_EXTRACTION_H
#define SPA_EXTRACTION_H

#include <src/metadata_table.h>
#include <src/image.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/single_particle/parallel_ft.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <vector>


class SpaExtraction
{
	public:
		
		// @TODO: replace output type by vector<BufferedImage>
		template <typename T>
		static std::vector<std::vector<Image<Complex>>> extractMovieStackFS(
				const MetaDataTable& mdt, 
				const RawImage<T>& movie,
				int boxSize,
				double outPs, double coordsPs, double moviePs, double dataPs, // @TODO: replace by data type
				const std::vector<std::vector<gravis::d2Vector>>* offsets_in,
				std::vector<std::vector<gravis::d2Vector>>* offsets_out,
				int num_threads);
};

template <typename T>
std::vector<std::vector<Image<Complex>>> SpaExtraction::extractMovieStackFS(
		const MetaDataTable& mdt, 
		const RawImage<T>& movie,
		int boxSize,
		double outPs, double coordsPs, double moviePs, double dataPs,
		const std::vector<std::vector<gravis::d2Vector>>* offsets_in,
		std::vector<std::vector<gravis::d2Vector>>* offsets_out,
		int num_threads)
{
	std::vector<std::vector<Image<Complex>>> out(mdt.numberOfObjects());
	const long pc = mdt.numberOfObjects();

	const int w0 = movie.xdim;
	const int h0 = movie.ydim;
	const int fc = movie.zdim;

	if (dataPs < 0) 
	{
		dataPs = outPs;
	}

	for (long p = 0; p < pc; p++)
	{
		out[p] = std::vector<Image<Complex>>(fc);
	}

	const int sqMg = 2*(int)(0.5 * boxSize * outPs / moviePs + 0.5);

	std::vector<ParFourierTransformer> fts(num_threads);

	std::vector<Image<RFLOAT>> aux0(num_threads);
	std::vector<Image<Complex>> aux1(num_threads);

	for (int th = 0; th < num_threads; th++)
	{
		aux0[th] = Image<RFLOAT>(sqMg, sqMg);

		if (outPs != moviePs)
		{
			aux1[th] = Image<Complex>(sqMg/2 + 1, sqMg);
		}
	}

	#pragma omp parallel for num_threads(num_threads)
	for (long int f = 0; f < fc; f++)
	{
		int tf = omp_get_thread_num();

		for (long p = 0; p < pc; p++)
		{
			int t = tf;

			out[p][f] = Image<Complex>(sqMg,sqMg);

			double xpC, ypC;

			mdt.getValue(EMDL_IMAGE_COORD_X, xpC, p);
			mdt.getValue(EMDL_IMAGE_COORD_Y, ypC, p);

			const double xpO = (int)(coordsPs * xpC / dataPs);
			const double ypO = (int)(coordsPs * ypC / dataPs);

			int x0 = (int)round(xpO * dataPs / moviePs) - sqMg / 2;
			int y0 = (int)round(ypO * dataPs / moviePs) - sqMg / 2;

			if (offsets_in != 0 && offsets_out != 0)
			{
				double dxM = (*offsets_in)[p][f].x * outPs / moviePs;
				double dyM = (*offsets_in)[p][f].y * outPs / moviePs;

				int dxI = (int)round(dxM);
				int dyI = (int)round(dyM);

				x0 += dxI;
				y0 += dyI;

				double dxR = (dxM - dxI) * moviePs / outPs;
				double dyR = (dyM - dyI) * moviePs / outPs;

				(*offsets_out)[p][f] = gravis::d2Vector(dxR, dyR);
			}

			for (long int y = 0; y < sqMg; y++)
			for (long int x = 0; x < sqMg; x++)
			{
				int xx = x0 + x;
				int yy = y0 + y;

				if (xx < 0) xx = 0;
				else if (xx >= w0) xx = w0 - 1;

				if (yy < 0) yy = 0;
				else if (yy >= h0) yy = h0 - 1;

				DIRECT_NZYX_ELEM(aux0[t].data, 0, 0, y, x) = movie(xx,yy,f);
			}

			if (outPs == moviePs)
			{
				fts[t].FourierTransform(aux0[t](), out[p][f]());
			}
			else
			{
				fts[t].FourierTransform(aux0[t](), aux1[t]());
				out[p][f] = FilterHelper::cropCorner2D(aux1[t], boxSize/2+1, boxSize);
			}

			out[p][f](0,0) = Complex(0.0,0.0);
		}
	}

	return out;
}


#endif
