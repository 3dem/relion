#ifndef JAZ_MOVIE_LOADER_H
#define JAZ_MOVIE_LOADER_H

#include <src/image.h>
#include <src/jaz/image/buffered_image.h>
#include <src/renderEER.h>
#include <string>

class MovieLoader
{
	public:

		template <typename T>
		static BufferedImage<T> readDense(
				std::string movieFn,
				const RawImage<RFLOAT>* gainRef,
				const RawImage<bool>* defectivePixels,
				int frame0,
				int numFrames,
				RFLOAT hot,
				int num_threads);

		template <typename T>
		static BufferedImage<T> readEER(
				std::string movieFn,
				const RawImage<RFLOAT>* gainRef,
				const RawImage<bool>* defectivePixels,
				int frame0,
				int numFrames,
				int eer_upsampling,
				int eer_grouping,
				int num_threads);

		template <typename T>
		static void fixDefects(
				RawImage<T>& muGraphFrame,
				const RawImage<bool>* defectivePixels,
				int num_threads, bool isEER);
};

template <typename T>
BufferedImage<T> MovieLoader::readDense(
			std::string movieFn,
			const RawImage<RFLOAT>* gainRef,
			const RawImage<bool>* defectivePixels,
			int frame0,
			int numFrames,
			RFLOAT hot,
			int num_threads)
{
	Image<float> mgStack;
	mgStack.read(movieFn, false, -1, false, true); // final true means 2D movies, not 3D map

	const std::string tag = "MovieLoader::readDense: ";

	if (mgStack.data.zdim > 1)
	{
		REPORT_ERROR_STR(tag << "the file " << movieFn << " looks like a 3D image and not a stack "
						 << "(the slices are along the Z dimension instead of N)");
	}

	const long int w0 = mgStack.data.xdim;
	const long int h0 = mgStack.data.ydim;
	const long int pixCt = w0 * h0;
	const int fc = numFrames;


	if (mgStack.data.ndim < frame0 + numFrames)
	{
		REPORT_ERROR_STR(
			tag << "insufficient number of frames in " << movieFn
			<< " (found: " << mgStack.data.ndim << ", expected: " << frame0 + numFrames << ")");
	}

	const bool useGain = gainRef != 0;
	if (useGain && (w0 != gainRef->xdim || h0 != gainRef->ydim))
	{
		REPORT_ERROR_STR(tag << "incompatible gain reference - size is different from " << movieFn);
	}

	const bool do_fixDefect = defectivePixels != 0;
	if (do_fixDefect && (w0 != defectivePixels->xdim || h0 != defectivePixels->ydim))
	{
		REPORT_ERROR_STR(tag << "incompatible defect mask - size is different from " << movieFn);
	}


	BufferedImage<T> out(w0, h0, fc);

	for (long f = 0; f < fc; f++)
	{
		Image<float> muGraphFrame_xmipp;
		muGraphFrame_xmipp.read(movieFn, true, frame0 + f, false, true);

		RawImage<T> muGraphFrame(muGraphFrame_xmipp);

		#pragma omp parallel for num_threads(num_threads)
		for (long int i = 0; i < pixCt; i++)
		{
			RFLOAT val = muGraphFrame[i];
			const RFLOAT gain = useGain? (*gainRef)[i] : 1;

			if (hot > 0.0 && val > hot) val = hot;

			 muGraphFrame[i] = -gain * val;
		}

		if (do_fixDefect)
		{
			fixDefects(muGraphFrame, defectivePixels, num_threads, false);
		}

		out.getSliceRef(f).copyFrom(muGraphFrame);
	}

	return out;
}

template <typename T>
BufferedImage<T> MovieLoader::readEER(
		std::string movieFn,
		const RawImage<RFLOAT>* gainRef,
		const RawImage<bool>* defectivePixels,
		int frame0,
		int numFrames,
		int eer_upsampling,
		int eer_grouping,
		int num_threads)
{
	EERRenderer renderer;
	renderer.read(movieFn, eer_upsampling);

	const long int w0 = renderer.getWidth();
	const long int h0 = renderer.getHeight();
	const long int pixCt = w0 * h0;
	const int fc = numFrames;
	const std::string tag = "MovieLoader::readDense: ";


	const bool useGain = gainRef != 0;
	if (useGain && (w0 != gainRef->xdim || h0 != gainRef->ydim))
	{
		REPORT_ERROR_STR(tag << "incompatible gain reference - size (x = " <<  gainRef->xdim <<
                                 ", y = " << gainRef->ydim << ") is different from " << movieFn <<
                                 " (x = " << w0 << ", y = " << h0 << ")");
	}

	const bool do_fixDefect = defectivePixels != 0;
	if (do_fixDefect && (w0 != defectivePixels->xdim || h0 != defectivePixels->ydim))
	{
		REPORT_ERROR_STR(tag << "incompatible defect mask - size (x = " << defectivePixels->xdim <<
                                 ", y = " << defectivePixels->ydim << ") is different from " << movieFn <<
                                 " (x = " << w0 << ", y = " << h0 << ")");
	}


	BufferedImage<T> out(w0, h0, fc);

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		MultidimArray<float> muGraphFrame_xmipp;
		// this takes 1-indexed frame numbers
		renderer.renderFrames(
					(frame0 + f) * eer_grouping + 1,
					(frame0 + f + 1) * eer_grouping,
					muGraphFrame_xmipp);

		RawImage<T> muGraphFrame(muGraphFrame_xmipp);

		if (useGain)
		{
			#pragma omp parallel for num_threads(num_threads)
			for (long int i = 0; i < pixCt; i++)
			{
				const RFLOAT val = muGraphFrame[i];
				const RFLOAT gain = useGain? (*gainRef)[i] : 1.0;

				muGraphFrame[i] = -gain * val;
			}
		}

		if (do_fixDefect)
		{
			fixDefects(muGraphFrame, defectivePixels, num_threads, true);
		}

		out.getSliceRef(f).copyFrom(muGraphFrame);
	}

	return out;
}

template <typename T>
void MovieLoader::fixDefects(
		RawImage<T>& muGraphFrame,
		const RawImage<bool>* defectivePixels,
		int num_threads, bool isEER)
{
	const long int w0 = muGraphFrame.xdim;
	const long int h0 = muGraphFrame.ydim;
	const long int pixCt = w0 * h0;

	RFLOAT frame_mean = 0;
	long int n_valid = 0;

	#pragma omp parallel for reduction(+:frame_mean, n_valid) num_threads(num_threads)
	for (long int n = 0; n < pixCt; ++n)
	{
		if (!(*defectivePixels)[n])
		{
			frame_mean += muGraphFrame[n];
			n_valid ++;
		}
	}

	frame_mean /= n_valid;

	RFLOAT frame_var = 0;

	#pragma omp parallel for reduction(+:frame_var) num_threads(num_threads)
	for (long int n = 0; n < pixCt; ++n)
	{
		if (!(*defectivePixels)[n])
		{
			RFLOAT d = (muGraphFrame[n] - frame_mean);
			frame_var += d * d;
		}
	}

	const RFLOAT frame_std = std::sqrt(frame_var / (n_valid - 1));

	const int min_num_ok = 6;
	const int d_max = isEER ? 4: 2;
	const int PBUF_SIZE = 100;

	#pragma omp parallel for num_threads(num_threads)
	for (long int y = 0; y < h0; y++)
	for (long int x = 0; x < w0; x++)
	{
		if (!(*defectivePixels)(x,y)) continue;

		int n_ok = 0;
		RFLOAT pbuf[PBUF_SIZE];

		for (int dy = -d_max; dy <= d_max; dy++)
		for (int dx = -d_max; dx <= d_max; dx++)
		{
			const int yy = y + dy;
			const int xx = x + dx;

			if (xx < 0 || xx >= w0 || yy < 0 || yy >= h0 || (*defectivePixels)(xx,yy))
				continue;

			pbuf[n_ok] = muGraphFrame(xx,yy);
			n_ok++;
		}
//						std::cout << "n_ok = " << n_ok << std::endl;
		if (n_ok > min_num_ok)
		{
			muGraphFrame(x,y) = pbuf[rand() % n_ok];
		}
		else
		{
			muGraphFrame(x,y) = rnd_gaus(frame_mean, frame_std);
		}
	}
}

#endif
