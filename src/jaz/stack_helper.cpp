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
#include <src/jaz/filter_helper.h>
#include <src/jaz/Fourier_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image_log.h>
#include <src/fftw.h>
#include <src/micrograph_model.h>
#include <src/jaz/resampling_helper.h>
#include <src/jaz/parallel_ft.h>

#include <omp.h>

using namespace gravis;

std::vector<MetaDataTable> StackHelper::splitByMicrographName(const MetaDataTable* mdt)
{
	std::vector<MetaDataTable> out(0);
	
	if (!mdt->labelExists(EMDL_MICROGRAPH_NAME))
	{
		REPORT_ERROR("StackHelper::splitByMicrographName: "+EMDL::label2Str(EMDL_MICROGRAPH_NAME)+" missing in meta_data_table.\n");
	}
	
	std::string testString;
	mdt->getValue(EMDL_MICROGRAPH_NAME, testString, 0);
	
	MetaDataTable md2(*mdt);
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
		const MetaDataTable *mdt, std::string path,
		int threads, std::vector<ParFourierTransformer> *fts,
		bool centerParticle)
{
	std::vector<Image<Complex> > out(mdt->numberOfObjects());
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
		int threadnum = omp_get_thread_num();
		
		std::string sliceName;
		mdt->getValue(EMDL_IMAGE_NAME, sliceName, i);
		Image<RFLOAT> in;
		in.read(sliceName, true, -1, false, true);
		
		(*fts)[threadnum].FourierTransform(in(), out[i]());
		
		if (centerParticle)
		{
			const int s = in.data.ydim;
			
			double xoff, yoff;
		
			mdt->getValue(EMDL_ORIENT_ORIGIN_X, xoff, i);
			mdt->getValue(EMDL_ORIENT_ORIGIN_Y, yoff, i);
			
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

std::vector<std::vector<Image<Complex>>> StackHelper::loadMovieStackFS(
		const MetaDataTable* mdt,
		std::string moviePath,
		bool center, int threads,
		std::vector<ParFourierTransformer> *fts,
		int firstFrame, int lastFrame, bool verbose)
{
	std::vector<std::vector<Image<Complex> > > out(mdt->numberOfObjects());
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
	
	if (verbose)
		std::cout << "loading " << finName << "\n";
	
	Image<RFLOAT> in;
	in.read(finName);
	
	if (verbose)
	{
		std::cout << "size = " << in.data.xdim << "x" << in.data.ydim << "x" << in.data.zdim << "x" << in.data.ndim << "\n";
		std::cout << "pc = " << pc << "\n";
	}
	
	const int fcMov = in.data.ndim / pc;
	
	if (fcMov < lastFrame+1)
	{
		std::stringstream sts;
		sts << "StackHelper::loadMovieStackFS: unable to load " << (lastFrame - firstFrame + 1) << " frames.";
		sts << "Only " << fcMov << " frames present in " << finName << ".\n";
		REPORT_ERROR(sts.str());
	}
	
	const int lfAct = lastFrame > 0? lastFrame : fcMov - 1;
	const int fc = lfAct - firstFrame + 1;
	
	const int w = in.data.xdim;
	const int h = in.data.ydim;
	
	for (long p = 0; p < pc; p++)
	{
		out[p] = std::vector<Image<Complex> >(fc);
		
		if (threads > 1)
		{
			#pragma omp parallel for num_threads(threads)
			for (long f = 0; f < fc; f++)
			{
				int threadnum = omp_get_thread_num();
				
				Image<RFLOAT> aux(w,h);
				
				SliceHelper::extractStackSlice(in, aux, (f+firstFrame)*pc + p);
				if (center) CenterFFT(aux(), true);
				(*fts)[threadnum].FourierTransform(aux(), out[p][f]());
			}
		}
		else
		{
			Image<RFLOAT> aux(w,h);
			FourierTransformer ft;
			
			for (long f = 0; f < fc; f++)
			{
				SliceHelper::extractStackSlice(in, aux, (f+firstFrame)*pc + p);
				if (center) CenterFFT(aux(), true);
				ft.FourierTransform(aux(), out[p][f]());
			}
		}
	}
	
	return out;
}

// deprecated:
std::vector<std::vector<Image<Complex>>> StackHelper::extractMovieStackFS(
		const MetaDataTable* mdt,
		std::string metaPath, std::string moviePath, std::string movie_ending,
		double outPs, double coordsPs, double moviePs,
		int squareSize, int threads,
		bool loadData, int firstFrame, int lastFrame,
		RFLOAT hot, bool verbose)
{
	std::vector<std::vector<Image<Complex>>> out(mdt->numberOfObjects());
	const long pc = mdt->numberOfObjects();
	
	std::string name0;
	mdt->getValue(EMDL_MICROGRAPH_NAME, name0, 0);
	
	std::string ending;
	
	if (movie_ending == "")
	{
		ending = name0.substr(name0.find_last_of(".")+1);
		
		if (verbose)
		{
			std::cout << "orig. name: " << name0 << "\n";
			std::cout << "ending: " << ending << "\n";
		}
	}
	else
	{
		if (verbose)
		{
			std::cout << "orig. name: " << name0 << "\n";
		}
		
		ending = movie_ending;
		name0 = name0.substr(0, name0.find_last_of(".")+1) + movie_ending;
		
		if (verbose)
		{
			std::cout << "act. name:  " << name0 << "\n";
		}
	}
	
	Image<float> mgStack;
	Image<RFLOAT> gainRef;
	
	std::string mgNameAct;
	
	bool useGain = false;
	
	if (ending == "star")
	{
		std::string metaNameAct;
		
		if (metaPath == "")
		{
			metaNameAct = name0;
		}
		else
		{
			metaNameAct = metaPath + "/" + name0.substr(name0.find_last_of("/")+1);
		}
		
		if (verbose)
		{
			std::cout << "loading: " << metaNameAct << "\n";
		}
		
		MetaDataTable metaMdt;
		metaMdt.read(metaNameAct, "general");
		
		std::string mgName0;
		
		if (!metaMdt.containsLabel(EMDL_MICROGRAPH_MOVIE_NAME))
		{
			REPORT_ERROR("StackHelper::extractMovieStackFS: " + metaNameAct
						 + " does not contain a rlnMicrographMovieName field.");
		}
		
		metaMdt.getValue(EMDL_MICROGRAPH_MOVIE_NAME, mgName0, 0);
		
		if (movie_ending != "")
		{
			mgName0.substr(0, mgName0.find_last_of(".")+1) + movie_ending;
		}
		
		std::string gainNameAct;
		
		if (moviePath == "")
		{
			mgNameAct = mgName0;
		}
		else
		{
			mgNameAct = moviePath + "/" + mgName0.substr(mgName0.find_last_of("/")+1);
		}
		
		if (verbose)
		{
			std::cout << "loading: " << mgNameAct << "\n";
		}
		
		mgStack.read(mgNameAct, false);
		
		if (loadData)
		{
			if (!metaMdt.containsLabel(EMDL_MICROGRAPH_GAIN_NAME))
			{
				REPORT_ERROR("StackHelper::extractMovieStackFS: " + metaNameAct
							 + " does not contain a rlnMicrographGainName field.");
			}
			
			std::string gainName0;
			metaMdt.getValue(EMDL_MICROGRAPH_GAIN_NAME, gainName0, 0);
			
			if (moviePath == "")
			{
				gainNameAct = gainName0;
			}
			else
			{
				gainNameAct = moviePath + "/" + gainName0.substr(gainName0.find_last_of("/")+1);
			}
			
			gainRef.read(gainNameAct, true);
			useGain = true;
			
			if (gainRef.data.xdim != mgStack.data.xdim
					|| gainRef.data.ydim != mgStack.data.ydim)
			{
				std::stringstream stsm;
				stsm << mgStack.data.xdim << "x" << mgStack.data.ydim;
				std::stringstream stsg;
				stsg << gainRef.data.xdim << "x" << gainRef.data.ydim;
				
				REPORT_ERROR("StackHelper::extractMovieStackFS: gain reference ("
							 + gainNameAct + ") and movie stack (" + mgNameAct
							 + ") are of unequal size ("+stsg.str()+" vs. "+stsm.str()+").");
			}
		}
	}
	else
	{
		if (moviePath == "")
		{
			mgNameAct = name0;
		}
		else
		{
			mgNameAct = moviePath + "/" + name0.substr(name0.find_last_of("/")+1);
		}
		
		if (verbose)
		{
			std::cout << name0 << " -> " << mgNameAct << "\n";
			std::cout << "loading: " << mgNameAct << "\n";
		}
		
		mgStack.read(mgNameAct, false);
	}
	
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
	const int fc = lastFrame > 0? lastFrame - firstFrame + 1 : fcM - firstFrame;
	
	if (verbose)
	{
		if (dataInZ) std::cout << "data in Z\n";
		else std::cout << "data in N\n";
		
		std::cout << "frame count in movie = " << fcM << "\n";
		std::cout << "frame count to load  = " << fc << "\n";
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
	
	Image<float> muGraph;
	
	for (long f = 0; f < fc; f++)
	{
		muGraph.read(mgNameAct, true, f+firstFrame);
		
		if (verbose) std::cout << (f+1) << "/" << fc << "\n";
		
		#pragma omp parallel for num_threads(threads)
		for (long p = 0; p < pc; p++)
		{
			int t = omp_get_thread_num();
			
			out[p][f] = Image<Complex>(sqMg,sqMg);
			
			double xpC, ypC;
			mdt->getValue(EMDL_IMAGE_COORD_X, xpC, p);
			mdt->getValue(EMDL_IMAGE_COORD_Y, ypC, p);
			
			const double xpO = (int)(coordsPs * xpC / outPs) - squareSize/2;
			const double ypO = (int)(coordsPs * ypC / outPs) - squareSize/2;
			
			const int x0 = (int)round(xpO * outPs / moviePs);
			const int y0 = (int)round(ypO * outPs / moviePs);
			
			// @TODO: read whole-image shifts from micrograph class
			
			for (long int y = 0; y < sqMg; y++)
			for (long int x = 0; x < sqMg; x++)
			{
				int xx = x0 + x;
				int yy = y0 + y;
				
				if (xx < 0) xx = 0;
				else if (xx >= w0) xx = w0 - 1;
				
				if (yy < 0) yy = 0;
				else if (yy >= h0) yy = h0 - 1;
				
				RFLOAT val = DIRECT_NZYX_ELEM(muGraph.data, 0, 0, yy, xx);
				RFLOAT gain = 1.0;
				
				if (useGain) gain = DIRECT_NZYX_ELEM(gainRef.data, 0, 0, yy, xx);
				
				if (hot > 0.0 && val > hot) val = hot;
				
				DIRECT_NZYX_ELEM(aux0[t].data, 0, 0, y, x) = -gain * val;
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


std::vector<std::vector<Image<Complex>>> StackHelper::extractMovieStackFS(
		const MetaDataTable* mdt,
		Image<RFLOAT>* gainRef, std::string movieFn,
		double outPs, double coordsPs, double moviePs,
		int squareSize, int threads,
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
		for (long p = 0; p < pc; p++)
		{
			int tp = omp_get_thread_num();
			
			int t = saveMemory? tp : tf;
			
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
				
				RFLOAT val = DIRECT_NZYX_ELEM(muGraph.data, 0, 0, yy, xx);
				RFLOAT gain = 1.0;
				
				if (useGain) gain = DIRECT_NZYX_ELEM(gainRef->data, 0, 0, yy, xx);
				
				if (hot > 0.0 && val > hot) val = hot;
				
				DIRECT_NZYX_ELEM(aux0[t].data, 0, 0, y, x) = -gain * val;
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

Image<Complex> StackHelper::projectView(Projector* projector, const MetaDataTable* mdt, int index)
{
	const int s = projector->ori_size;
	const int sh = s/2 + 1;
	
	double xoff, yoff;
	
	mdt->getValue(EMDL_ORIENT_ORIGIN_X, xoff, index);
	mdt->getValue(EMDL_ORIENT_ORIGIN_Y, yoff, index);
	
	double rot, tilt, psi;
	
	Matrix2D<RFLOAT> A3D;
	mdt->getValue(EMDL_ORIENT_ROT, rot, index);
	mdt->getValue(EMDL_ORIENT_TILT, tilt, index);
	mdt->getValue(EMDL_ORIENT_PSI, psi, index);
	Euler_angles2matrix(rot, tilt, psi, A3D);
	
	Image<Complex> out(sh, s);
	out.data.initZeros();
	
	projector->get2DFourierTransform(out.data, A3D, false);
	shiftImageInFourierTransform(out(), out(), s, s/2 - xoff, s/2 - yoff);
	
	return out;
}

std::vector<Image<Complex> > StackHelper::projectStack(Projector* projector, const MetaDataTable* mdt)
{
	std::vector<Image<Complex> > out(mdt->numberOfObjects());
	const long ic = mdt->numberOfObjects();
	
	const int s = projector->ori_size;
	const int sh = s/2 + 1;
	
	for (long i = 0; i < ic; i++)
	{
		double xoff, yoff;
		
		mdt->getValue(EMDL_ORIENT_ORIGIN_X, xoff, i);
		mdt->getValue(EMDL_ORIENT_ORIGIN_Y, yoff, i);
		
		double rot, tilt, psi;
		
		Matrix2D<RFLOAT> A3D;
		mdt->getValue(EMDL_ORIENT_ROT, rot, i);
		mdt->getValue(EMDL_ORIENT_TILT, tilt, i);
		mdt->getValue(EMDL_ORIENT_PSI, psi, i);
		Euler_angles2matrix(rot, tilt, psi, A3D);
		
		out[i] = Image<Complex>(sh, s);
		out[i].data.initZeros();
		projector->get2DFourierTransform(out[i].data, A3D, false);
		shiftImageInFourierTransform(out[i](), out[i](), s, s/2 - xoff, s/2 - yoff);
	}
	
	return out;
}

std::vector<Image<Complex> > StackHelper::projectStackPar(
		Projector* projector, const MetaDataTable* mdt, int numThreads)
{
	std::vector<Image<Complex> > out(mdt->numberOfObjects());
	const long ic = mdt->numberOfObjects();
	
	const int s = projector->ori_size;
	const int sh = s/2 + 1;
	
	#pragma omp parallel for num_threads(numThreads)
	for (long i = 0; i < ic; i++)
	{
		double xoff, yoff;
		
		mdt->getValue(EMDL_ORIENT_ORIGIN_X, xoff, i);
		mdt->getValue(EMDL_ORIENT_ORIGIN_Y, yoff, i);
		
		double rot, tilt, psi;
		
		Matrix2D<RFLOAT> A3D;
		mdt->getValue(EMDL_ORIENT_ROT, rot, i);
		mdt->getValue(EMDL_ORIENT_TILT, tilt, i);
		mdt->getValue(EMDL_ORIENT_PSI, psi, i);
		Euler_angles2matrix(rot, tilt, psi, A3D);
		
		out[i] = Image<Complex>(sh, s);
		out[i].data.initZeros();
		projector->get2DFourierTransform(out[i].data, A3D, false);
		shiftImageInFourierTransform(out[i](), out[i](), s, s/2 - xoff, s/2 - yoff);
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
	
	std::cout << ic << " images in stack.\n";
	
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

std::vector<Image<Complex> > StackHelper::applyBeamTilt(std::vector<Image<Complex> > &stack,
														RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
														RFLOAT tilt_x, RFLOAT tilt_y,
														RFLOAT tilt_xx, RFLOAT tilt_xy, RFLOAT tilt_yy)
{
	std::vector<Image<Complex> > out(stack.size());
	const long ic = stack.size();
	
	for (long i = 0; i < ic; i++)
	{
		out[i] = stack[i];
		selfApplyBeamTilt(
					out[i].data, tilt_x, tilt_y,
					tilt_xx, tilt_xy, tilt_yy,
					lambda, Cs, angpix, stack[i].data.ydim);
	}
	
	return out;
	
}

std::vector<Image<Complex> > StackHelper::applyBeamTiltPar(std::vector<Image<Complex> > &stack,
														   RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
														   RFLOAT tilt_x, RFLOAT tilt_y,
														   RFLOAT tilt_xx, RFLOAT tilt_xy, RFLOAT tilt_yy, int numThreads)
{
	std::vector<Image<Complex> > out(stack.size());
	const long ic = stack.size();
	
	#pragma omp parallel for num_threads(numThreads)
	for (long i = 0; i < ic; i++)
	{
		out[i] = stack[i];
		selfApplyBeamTilt(
					out[i].data, tilt_x, tilt_y,
					tilt_xx, tilt_xy, tilt_yy,
					lambda, Cs, angpix, stack[i].data.ydim);
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
