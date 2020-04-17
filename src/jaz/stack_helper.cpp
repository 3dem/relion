/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <src/jaz/stack_helper.h>
#include <src/jaz/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/Fourier_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image_log.h>
#include <src/fftw.h>
#include <src/micrograph_model.h>
#include <src/jaz/resampling_helper.h>
#include <src/jaz/parallel_ft.h>
#include <src/jaz/new_ft.h>

#include <omp.h>

using namespace gravis;

std::vector<MetaDataTable> StackHelper::splitByMicrographName(const MetaDataTable& mdt)
{
	std::vector<MetaDataTable> out(0);

	if (!mdt.labelExists(EMDL_MICROGRAPH_NAME))
	{
		REPORT_ERROR("StackHelper::splitByMicrographName: "
					 + EMDL::label2Str(EMDL_MICROGRAPH_NAME)
					 + " missing from MetaDataTable.\n");
	}

	MetaDataTable md2(mdt);
	md2.newSort(EMDL_MICROGRAPH_NAME);

	const long lc = md2.numberOfObjects();
	std::string lastName = "", curName;
	long curInd = -1;

	for (int i = 0; i < lc; i++)
	{
		md2.getValue(EMDL_MICROGRAPH_NAME, curName, i);

		if (curName != lastName)
		{
			lastName = curName;
			curInd++;
			out.push_back(MetaDataTable());
		}

		out[curInd].addObject(md2.getObject(i));
	}

	for (int i = 0; i <= curInd; i++)
	{
		out[i].newSort(EMDL_IMAGE_NAME, false, false, true);
	}

	return out;
}

MetaDataTable StackHelper::merge(const std::vector<MetaDataTable> &mdts)
{
	MetaDataTable out;

	for (int i = 0; i < mdts.size(); i++)
	{
		out.append(mdts[i]);
	}

	return out;
}

std::vector<MetaDataTable> StackHelper::splitByStack(const MetaDataTable* mdt)
{
	std::vector<MetaDataTable> out(0);

	if (!mdt->labelExists(EMDL_IMAGE_NAME))
	{
		REPORT_ERROR("StackHelper::splitByStack: "+EMDL::label2Str(EMDL_IMAGE_NAME)+" missing in meta_data_table.\n");
	}

	std::string testString;
	mdt->getValue(EMDL_IMAGE_NAME, testString, 0);

	if (testString.find("@") < 0)
	{
		REPORT_ERROR("StackHelper::splitByStack: "+EMDL::label2Str(EMDL_IMAGE_NAME)+" does not contain an '@'.\n");
	}

	MetaDataTable md2(*mdt);
	md2.newSort(EMDL_IMAGE_NAME, false, true);

	const long lc = md2.numberOfObjects();
	std::string lastName = "", curName, curFullName;
	long curInd = -1;

	for (int i = 0; i < lc; i++)
	{
		md2.getValue(EMDL_IMAGE_NAME, curFullName, i);

		curName = curFullName.substr(curFullName.find("@")+1);

		if (curName != lastName)
		{
			lastName = curName;
			curInd++;
			out.push_back(MetaDataTable());
		}

		out[curInd].addObject(md2.getObject(i));
	}

	for (int i = 0; i <= curInd; i++)
	{
		out[i].newSort(EMDL_IMAGE_NAME, false, false, true);
	}

	return out;
}

std::vector<Image<RFLOAT> > StackHelper::loadStack(const MetaDataTable* mdt, std::string path, int threads)
{
	std::vector<Image<RFLOAT>> out(mdt->numberOfObjects());
	const long ic = mdt->numberOfObjects();

	std::string name, fullName;
	mdt->getValue(EMDL_IMAGE_NAME, fullName, 0);
	name = fullName.substr(fullName.find("@")+1);

	if (path != "")
	{
		name = path + "/" + name.substr(name.find_last_of("/")+1);
	}

	#pragma omp parallel for num_threads(threads)
	for (long i = 0; i < ic; i++)
	{
		std::string sliceName;
		mdt->getValue(EMDL_IMAGE_NAME, sliceName, i);	
		out[i].read(sliceName, true, -1, false, true);
	}

	return out;
}

std::vector<Image<Complex> > StackHelper::loadStackFS(
		const MetaDataTable& mdt, std::string path,
		int threads, bool centerParticle, ObservationModel* obs)
{
	std::vector<Image<Complex> > out(mdt.numberOfObjects());

	if (centerParticle && obs == 0)
	{
		REPORT_ERROR("StackHelper::loadStackFS: centering particles requires an observation model.");
	}

	const long ic = mdt.numberOfObjects();

	std::string name, fullName;
	mdt.getValue(EMDL_IMAGE_NAME, fullName, 0);
	name = fullName.substr(fullName.find("@")+1);

	if (path != "")
	{
		name = path + "/" + name.substr(name.find_last_of("/")+1);
	}

	Image<RFLOAT> dummy;
	dummy.read(name, false);

	const int s = dummy.data.xdim;

	NewFFTPlan<RFLOAT>::type plan(s,s,1);

	#pragma omp parallel for num_threads(threads)
	for (long i = 0; i < ic; i++)
	{
		int optGroup = obs->getOpticsGroup(mdt, i);
		double angpix = obs->getPixelSize(optGroup);

		std::string sliceName;
		mdt.getValue(EMDL_IMAGE_NAME, sliceName, i);
		Image<RFLOAT> in;
		in.read(sliceName, true, -1, false, true);
		
		NewFFT::FourierTransform(in(), out[i](), plan);

		if (centerParticle)
		{
			const int s = in.data.ydim;

			double xoff, yoff;

			mdt.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, i);
			mdt.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, i);

			xoff /= angpix;
			yoff /= angpix;

			shiftImageInFourierTransform(out[i](), out[i](), s, xoff - s/2, yoff - s/2);
		}
	}

	return out;
}

void StackHelper::saveStack(std::vector<Image<RFLOAT> > &stack, std::string fn)
{
	const int w = stack[0].data.xdim;
	const int h = stack[0].data.ydim;
	const int c = stack.size();

	Image<RFLOAT> img(w,h,1,c);

	for (int i = 0; i < c; i++)
	{
		SliceHelper::insertStackSlice(stack[i], img, i);
	}

	img.write(fn);
}

std::vector<std::vector<Image<RFLOAT> > > StackHelper::loadMovieStack(const MetaDataTable* mdt, std::string moviePath)
{
	std::vector<std::vector<Image<RFLOAT> > > out(mdt->numberOfObjects());
	const long pc = mdt->numberOfObjects();

	std::string name, fullName, movieName;
	mdt->getValue(EMDL_IMAGE_NAME, fullName, 0);
	mdt->getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
	name = fullName.substr(fullName.find("@")+1);

	std::string finName;

	if (moviePath == "")
	{
		finName = name;
	}
	else
	{
		finName = moviePath + "/" + movieName.substr(movieName.find_last_of("/")+1);
	}

	std::cout << "loading real: " << finName << "\n";

	Image<RFLOAT> in;
	in.read(finName);

	std::cout << "size = " << in.data.xdim << "x" << in.data.ydim << "x" << in.data.zdim << "x" << in.data.ndim << "\n";
	std::cout << "pc = " << pc << "\n";

	const int fc = in.data.ndim / pc;

	const int w = in.data.xdim;
	const int h = in.data.ydim;

	for (long p = 0; p < pc; p++)
	{
		out[p] = std::vector<Image<RFLOAT> >(fc);

		for (long f = 0; f < fc; f++)
		{
			out[p][f] = Image<RFLOAT>(w,h);
			SliceHelper::extractStackSlice(in, out[p][f], f*pc + p);
		}
	}

	return out;
}

std::vector<std::vector<Image<Complex>>> StackHelper::extractMovieStackFS(
		const MetaDataTable* mdt,
		Image<RFLOAT>* gainRef, MultidimArray<bool>* defectMask, std::string movieFn,
		double outPs, double coordsPs, double moviePs,
		int squareSize, int threads, // squareSize is the output box size in pixels after downsampling to outPs
		bool loadData, int firstFrame, int lastFrame,
		RFLOAT hot, bool verbose, bool saveMemory,
		const std::vector<std::vector<gravis::d2Vector>>* offsets_in,
		std::vector<std::vector<gravis::d2Vector>>* offsets_out)
{
	std::vector<std::vector<Image<Complex>>> out(mdt->numberOfObjects());
	const long pc = mdt->numberOfObjects();

	Image<float> mgStack;
	mgStack.read(movieFn, false);

	if (verbose)
	{
		std::cout << "size: "
				  << mgStack().xdim << "x"
				  << mgStack().ydim << "x"
				  << mgStack().zdim << "x"
				  << mgStack().ndim << "\n";
	}

	const bool dataInZ = mgStack.data.zdim > 1;

	const int w0 = mgStack.data.xdim;
	const int h0 = mgStack.data.ydim;
	const int fcM = dataInZ? mgStack.data.zdim : mgStack.data.ndim;
	// lastFrame and firstFrame is 0 indexed, while fcM is 1-indexed
	const int fc = lastFrame > 0? lastFrame - firstFrame + 1 : fcM - firstFrame;

	if (fcM <= lastFrame)
	{
		REPORT_ERROR("StackHelper::extractMovieStackFS: insufficient number of frames in "+movieFn);
	}

	const bool useGain = gainRef != 0;
	if (useGain && (w0 != gainRef->data.xdim || h0 != gainRef->data.ydim))
	{
		REPORT_ERROR("StackHelper::extractMovieStackFS: incompatible gain reference - size is different from "+movieFn);
	}

	const bool fixDefect = false; // TAKANORI DEBUG: defectMask != 0;
	if (fixDefect && (w0 != defectMask->xdim || h0 != defectMask->ydim))
	{
		REPORT_ERROR("StackHelper::extractMovieStackFS: incompatible defect mask - size is different from "+movieFn);
	}

	if (verbose)
	{
		if (dataInZ) std::cout << "data in Z\n";
		else std::cout << "data in N\n";

		std::cout << "frame count in movie = " << fcM << "\n";
		std::cout << "frame count to load  = " << fc << "\n";

		std::cout << "pc, fc = " << pc << ", " << fc << "\n";
	}

	for (long p = 0; p < pc; p++)
	{
		out[p] = std::vector<Image<Complex>>(fc);
	}

	if (!loadData) return out;

	const int sqMg = 2*(int)(0.5 * squareSize * outPs / moviePs + 0.5);

	if (verbose)
	{
		std::cout << "square size in micrograph: " << sqMg << "\n";
	}

	std::vector<ParFourierTransformer> fts(threads);

	std::vector<Image<RFLOAT>> aux0(threads);
	std::vector<Image<Complex>> aux1(threads);

	for (int t = 0; t < threads; t++)
	{
		aux0[t] = Image<RFLOAT>(sqMg, sqMg);

		if (outPs != moviePs)
		{
			aux1[t] = Image<Complex>(sqMg/2+1,sqMg);
		}
	}


	int threads_f = saveMemory? 1 : threads;
	int threads_p = saveMemory? threads : 1;

	#pragma omp parallel for num_threads(threads_f)
	for (long f = 0; f < fc; f++)
	{
		int tf = omp_get_thread_num();

		Image<float> muGraph;
		muGraph.read(movieFn, true, f+firstFrame, false, true);

		if (verbose) std::cout << (f+1) << "/" << fc << "\n";

		#pragma omp parallel for num_threads(threads_p)
		for (long int y = 0; y < h0; y++)
		for (long int x = 0; x < w0; x++)
		{
			RFLOAT val = DIRECT_NZYX_ELEM(muGraph.data, 0, 0, y, x);
			RFLOAT gain = 1.0;

			if (useGain) gain = DIRECT_NZYX_ELEM(gainRef->data, 0, 0, y, x);
			if (hot > 0.0 && val > hot) val = hot;

			 DIRECT_NZYX_ELEM(muGraph.data, 0, 0, y, x) = -gain * val;
		}

		if (fixDefect)
		{
			RFLOAT frame_mean = 0, frame_std = 0;
			long long n_valid = 0;

			#pragma omp parallel for reduction(+:frame_mean, n_valid) num_threads(threads_p)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(muGraph.data) {
				if (!DIRECT_MULTIDIM_ELEM(*defectMask, n)) continue;
				frame_mean += DIRECT_MULTIDIM_ELEM(muGraph.data, n);
				n_valid ++;
			}
			frame_mean /=  n_valid;

			#pragma omp parallel for reduction(+:frame_std) num_threads(threads_p)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(muGraph.data) {
				if (!DIRECT_MULTIDIM_ELEM(*defectMask, n)) continue;
				RFLOAT d = (DIRECT_MULTIDIM_ELEM(muGraph.data, n) - frame_mean);
				frame_std += d * d;
			}
			frame_std = std::sqrt(frame_std / n_valid);

			// 25 neighbours; should be enough even for super-resolution images.
			const int NUM_MIN_OK = 6;
			const int D_MAX = 2;
			#pragma omp parallel for num_threads(threads_p)
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(muGraph.data)
			{
				if (!DIRECT_A2D_ELEM(*defectMask, i, j)) continue;

				int n_ok = 0;
				RFLOAT val = 0;
				for (int dy= -D_MAX; dy <= D_MAX; dy++)
				{
					int y = i + dy;
					if (y < 0 || y >= h0) continue;
					for (int dx = -D_MAX; dx <= D_MAX; dx++)
					{
						int x = j + dx;
						if (x < 0 || x >= w0) continue;
						if (DIRECT_A2D_ELEM(*defectMask, y, x)) continue;

						n_ok++;
						val += DIRECT_A2D_ELEM(muGraph.data, y, x);
					}
				}
				if (n_ok > NUM_MIN_OK) DIRECT_A2D_ELEM(muGraph.data, i, j) = val / n_ok;
				else DIRECT_A2D_ELEM(muGraph.data, i, j) = rnd_gaus(frame_mean, frame_std);
			}
		}

		// TODO: TAKANORI: Cache muGraph HERE

		#pragma omp parallel for num_threads(threads_p)
		for (long p = 0; p < pc; p++)
		{
			int tp = omp_get_thread_num();

			int t = saveMemory? tp : tf;

			out[p][f] = Image<Complex>(sqMg,sqMg);

			double xpC, ypC;

			mdt->getValue(EMDL_IMAGE_COORD_X, xpC, p);
			mdt->getValue(EMDL_IMAGE_COORD_Y, ypC, p);

			// TAKANORI: TODO: BUG? Probably this causes 1 px error when outPS changes
			const double xpO = (int)(coordsPs * xpC / outPs) - squareSize/2;
			const double ypO = (int)(coordsPs * ypC / outPs) - squareSize/2;

			int x0 = (int)round(xpO * outPs / moviePs);
			int y0 = (int)round(ypO * outPs / moviePs);

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

				(*offsets_out)[p][f] = d2Vector(dxR, dyR);
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

				DIRECT_NZYX_ELEM(aux0[t].data, 0, 0, y, x) = DIRECT_NZYX_ELEM(muGraph.data, 0, 0, yy, xx);
			}

			if (outPs == moviePs)
			{
				fts[t].FourierTransform(aux0[t](), out[p][f]());
			}
			else
			{
				fts[t].FourierTransform(aux0[t](), aux1[t]());
				out[p][f] = FilterHelper::cropCorner2D(aux1[t], squareSize/2+1, squareSize);
			}

			out[p][f](0,0) = Complex(0.0,0.0);
		}
	}

	return out;
}

// TAKANORI: TODO: Code duplication with above will be sorted out later!
std::vector<std::vector<Image<Complex>>> StackHelper::extractMovieStackFS(
		const MetaDataTable* mdt, std::vector<MultidimArray<float> > &Iframes,
		double outPs, double coordsPs, double moviePs,
		int squareSize, int threads,
		bool loadData,
		bool verbose,
		const std::vector<std::vector<gravis::d2Vector>>* offsets_in,
		std::vector<std::vector<gravis::d2Vector>>* offsets_out)
{
	std::vector<std::vector<Image<Complex>>> out(mdt->numberOfObjects());
	const long pc = mdt->numberOfObjects();

	const int fc = Iframes.size();
	if (fc == 0)
		REPORT_ERROR("Empty Iframes passed to StackHelper::extractMovieStackFS");
	const int w0 = Iframes[0].xdim;
	const int h0 = Iframes[0].ydim;

	if (verbose)
	{
		std::cout << "pc, fc = " << pc << ", " << fc << "\n";
		std::cout << "size: x = " << w0 << " y = " << h0 << "\n";
	}

	for (long p = 0; p < pc; p++)
	{
		out[p] = std::vector<Image<Complex>>(fc);
	}

	if (!loadData) return out;

	const int sqMg = 2*(int)(0.5 * squareSize * outPs / moviePs + 0.5);

	if (verbose)
	{
		std::cout << "square size in micrograph: " << sqMg << "\n";
	}

	std::vector<ParFourierTransformer> fts(threads);

	std::vector<Image<RFLOAT>> aux0(threads);
	std::vector<Image<Complex>> aux1(threads);

	for (int t = 0; t < threads; t++)
	{
		aux0[t] = Image<RFLOAT>(sqMg, sqMg);

		if (outPs != moviePs)
		{
			aux1[t] = Image<Complex>(sqMg/2+1,sqMg);
		}
	}

	#pragma omp parallel for num_threads(threads)
	for (long f = 0; f < fc; f++)
	{
		int tf = omp_get_thread_num();

		if (verbose) std::cout << (f+1) << "/" << fc << "\n";

		for (long p = 0; p < pc; p++)
		{
			int t = tf;

			out[p][f] = Image<Complex>(sqMg,sqMg);

			double xpC, ypC;

			mdt->getValue(EMDL_IMAGE_COORD_X, xpC, p);
			mdt->getValue(EMDL_IMAGE_COORD_Y, ypC, p);

			const double xpO = (int)(coordsPs * xpC / outPs) - squareSize/2;
			const double ypO = (int)(coordsPs * ypC / outPs) - squareSize/2;

			int x0 = (int)round(xpO * outPs / moviePs);
			int y0 = (int)round(ypO * outPs / moviePs);

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

				(*offsets_out)[p][f] = d2Vector(dxR, dyR);
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

				// Note the MINUS here!!!
				DIRECT_NZYX_ELEM(aux0[t].data, 0, 0, y, x) = -DIRECT_A2D_ELEM(Iframes[f], yy, xx);
			}

			if (outPs == moviePs)
			{
				fts[t].FourierTransform(aux0[t](), out[p][f]());
			}
			else
			{
				fts[t].FourierTransform(aux0[t](), aux1[t]());
				out[p][f] = FilterHelper::cropCorner2D(aux1[t], squareSize/2+1, squareSize);
			}

			out[p][f](0,0) = Complex(0.0,0.0);
		}
	}

	return out;
}
std::vector<Image<Complex> > StackHelper::FourierTransform(std::vector<Image<RFLOAT> >& stack)
{
	std::vector<Image<Complex> > out(stack.size());
	const long ic = stack.size();

	for (long i = 0; i < ic; i++)
	{
		FourierTransformer ft;
		ft.FourierTransform(stack[i].data, out[i].data);
	}

	return out;
}

std::vector<Image<RFLOAT> > StackHelper::inverseFourierTransform(std::vector<Image<Complex> >& stack)
{
	std::vector<Image<RFLOAT> > out(stack.size());
	const long ic = stack.size();

	const int h = stack[0].data.ydim;
	const int ww = stack[0].data.xdim;
	const int w = 2*(ww - 1);

	for (long i = 0; i < ic; i++)
	{
		out[i] = Image<RFLOAT>(w,h);

		FourierTransformer ft;
		ft.inverseFourierTransform(stack[i].data, out[i].data);
	}

	return out;
}

Image<RFLOAT> StackHelper::toSingleImage(const std::vector<Image<RFLOAT>> stack)
{
	const int s = stack.size();

	if (s < 1) return Image<RFLOAT>(0,0,0);

	const int w = stack[0].data.xdim;
	const int h = stack[0].data.ydim;

	Image<RFLOAT> out(w,h,1,s);

	for (int n = 0; n < s; n++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			DIRECT_NZYX_ELEM(out(),n,0,y,x) = stack[n](y,x);
		}
	}

	return out;
}

void StackHelper::varianceNormalize(std::vector<Image<Complex>>& movie, bool circleCropped)
{
	const int fc = movie.size();
	const int w = movie[0].data.xdim;
	const int h = movie[0].data.ydim;
	const int wt = 2*(w-1);

	double var = 0.0;
	double cnt = 0.0;

	const double rr = (w-2)*(w-2);

	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			if (x == 0 && y == 0) continue;

			if (circleCropped)
			{
				const double yy = y < w? y : y - h;
				const double xx = x;

				if (xx*xx + yy*yy > rr) continue;
			}

			double scale = x > 0? 2.0 : 1.0;

			var += scale * movie[f](y,x).norm();
			cnt += scale;
		}
	}

	const double scale = sqrt(wt*h*var/(cnt*fc));

	//std::cout << "scale: " << scale << "\n";

	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
			for (int x = 0; x < w; x++)
			{
				movie[f](y,x) /= scale;
			}
	}
}

std::vector<double> StackHelper::powerSpectrum(const std::vector<std::vector<Image<Complex>>> &stack)
{
	const int ic = stack.size();
	const int fc = stack[0].size();
	const int w = stack[0][0].data.xdim;
	const int h = stack[0][0].data.ydim;

	std::vector<double> out(w, 0.0), wgh(w, 0.0);

	for (int i = 0; i < ic; i++)
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			const Complex z = DIRECT_A2D_ELEM(stack[i][f].data, y, x);

			const double yy = y < w? y : y - h;
			const double xx = x;

			const int r = (int) sqrt(xx*xx + yy*yy);

			if (r >= w) continue;

			out[r] += z.norm();
			wgh[r] += 1.0;
		}
	}

	for (int x = 0; x < w; x++)
	{
		if (wgh[x] > 0.0)
		{
			out[x] /= wgh[x];
		}
	}

	return out;
}

std::vector<double> StackHelper::varSpectrum(const std::vector<std::vector<Image<Complex>>> &stack)
{
	const int ic = stack.size();
	const int fc = stack[0].size();
	const int w = stack[0][0].data.xdim;
	const int h = stack[0][0].data.ydim;

	std::vector<double> out(w, 0.0), wgh(w, 0.0);
	std::vector<Complex> mean(w, 0.0);

	for (int i = 0; i < ic; i++)
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			const Complex z = DIRECT_A2D_ELEM(stack[i][f].data, y, x);

			const double yy = y < w? y : y - h;
			const double xx = x;

			const int r = (int) sqrt(xx*xx + yy*yy);

			if (r >= w) continue;

			mean[r] += z;
			wgh[r] += 1.0;
		}
	}

	for (int x = 0; x < w; x++)
	{
		if (wgh[x] > 0.0)
		{
			mean[x] /= wgh[x];
		}
	}

	for (int i = 0; i < ic; i++)
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			const Complex z = DIRECT_A2D_ELEM(stack[i][f].data, y, x);

			const double yy = y < w? y : y - h;
			const double xx = x;

			const int r = (int) sqrt(xx*xx + yy*yy);

			if (r >= w) continue;

			out[r] += (z - mean[r]).norm();
		}
	}

	for (int x = 0; x < w; x++)
	{
		if (wgh[x] > 1.0)
		{
			out[x] /= (wgh[x] - 1.0);
		}
	}

	return out;
}

std::vector<double> StackHelper::powerSpectrum(
		const std::vector<std::vector<Image<Complex>>>& obs,
		const std::vector<Image<Complex> >& signal)
{
	const int ic = obs.size();
	const int fc = obs[0].size();
	const int w = obs[0][0].data.xdim;
	const int h = obs[0][0].data.ydim;

	std::vector<double> out(w, 0.0), wgh(w, 0.0);

	for (int i = 0; i < ic; i++)
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			const Complex z = DIRECT_A2D_ELEM(obs[i][f].data, y, x) - DIRECT_A2D_ELEM(signal[i].data, y, x);

			const double yy = y < w? y : y - h;
			const double xx = x;

			const int r = (int) sqrt(xx*xx + yy*yy);

			if (r >= w) continue;

			out[r] += z.norm();
			wgh[r] += 1.0;
		}
	}

	for (int x = 0; x < w; x++)
	{
		if (wgh[x] > 0.0)
		{
			out[x] /= wgh[x];
		}
	}

	return out;
}
