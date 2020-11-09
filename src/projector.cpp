/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
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
#include "src/projector.h"
#include <src/jaz/gravis/t3Vector.h>
#include <src/time.h>
#include <src/jaz/image/buffered_image.h>
#ifdef CUDA
#include <cufft.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#define SEMKEY 826976737978L /* key value for semget() */
#define PERMS 0666
#endif
//#define DEBUG

//#define PROJ_TIMING
#ifdef PROJ_TIMING
	Timer proj_timer;
	int TIMING_TOP = proj_timer.setNew("PROJECTOR - computeFourierTransformMap");
	int TIMING_GRID = proj_timer.setNew("PROJECTOR - gridCorr");
	int TIMING_PAD = proj_timer.setNew("PROJECTOR - padTransMap");
	int TIMING_CENTER = proj_timer.setNew("PROJECTOR - centerFFT");
	int TIMING_TRANS = proj_timer.setNew("PROJECTOR - transform");
	int TIMING_FAUX = proj_timer.setNew("PROJECTOR - Faux");
	int TIMING_POW = proj_timer.setNew("PROJECTOR - power_spectrum");
	int TIMING_INIT1 = proj_timer.setNew("PROJECTOR - init1");
	int TIMING_INIT2 = proj_timer.setNew("PROJECTOR - init2");
#define TIMING_TIC(id) proj_timer.tic(id)
#define TIMING_TOC(id) proj_timer.toc(id)
#else
#define TIMING_TIC(id)
#define TIMING_TOC(id)
#endif

using namespace gravis;

void Projector::initialiseData(int current_size)
{
	// By default r_max is half ori_size
	if (current_size < 0)
		r_max = ori_size / 2;
	else
		r_max = current_size / 2;

	// Never allow r_max beyond Nyquist...
	r_max = XMIPP_MIN(r_max, ori_size / 2);

	// Set pad_size
	pad_size = 2 * (ROUND(padding_factor * r_max) + 1) + 1;

	// Short side of data array
	switch (ref_dim)
	{
	case 2:
		data.resize(pad_size, pad_size / 2 + 1);
		break;
	case 3:
		data.resize(pad_size, pad_size, pad_size / 2 + 1);
		break;
	default:
		REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");
	}

	// Set origin in the y.z-center, but on the left side for x.
	data.setXmippOrigin();
	data.xinit=0;

}
void Projector::initZeros(int current_size)
{
	initialiseData(current_size);
	data.initZeros();
}

long int Projector::getSize()
{
	// Short side of data array
	switch (ref_dim)
	{
		case 2:
			return pad_size * (pad_size / 2 + 1);
			break;
		case 3:
			return pad_size * pad_size * (pad_size / 2 + 1);
			break;
		default:
			REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");
	}

}


// Fill data array with oversampled Fourier transform, and calculate its power spectrum
void Projector::computeFourierTransformMap(
		MultidimArray<RFLOAT> &vol_in, MultidimArray<RFLOAT> &power_spectrum,
		int current_size, int nr_threads, bool do_gridding, bool do_heavy, int min_ires,
		const MultidimArray<RFLOAT>* fourier_mask, bool do_gpu)
{
	TIMING_TIC(TIMING_TOP);

	TIMING_TIC(TIMING_INIT1);
	MultidimArray<RFLOAT> Mpad;
	MultidimArray<Complex > Faux;
	FourierTransformer transformer;
	RFLOAT normfft;

	// Size of padded real-space volume
	int padoridim = ROUND(padding_factor * ori_size);
	// make sure padoridim is even
	padoridim += padoridim%2;
	// Re-calculate padding factor
	padding_factor = (float)padoridim/(float)ori_size;

	// Initialize data array of the oversampled transform
	ref_dim = vol_in.getDim();

	bool do_fourier_mask = (fourier_mask != NULL);

	// Make Mpad
	switch (ref_dim)
	{
	case 2:
		if(do_heavy)
			Mpad.initZeros(padoridim, padoridim);
		else
			Mpad.reshape(padoridim, padoridim);
		if (data_dim == 2)
			normfft = (RFLOAT)(padding_factor * padding_factor);
		else
			normfft = (RFLOAT)(padding_factor * padding_factor * ori_size);
		break;
	case 3:
		if(do_heavy)
			Mpad.initZeros(padoridim, padoridim, padoridim);
		else
			Mpad.reshape(padoridim, padoridim, padoridim);
		if (data_dim == 3)
			normfft = (RFLOAT)(padding_factor * padding_factor * padding_factor);
		else
			normfft = (RFLOAT)(padding_factor * padding_factor * padding_factor * ori_size);
		break;
	default:
		REPORT_ERROR("Projector::computeFourierTransformMap%%ERROR: Dimension of the data array should be 2 or 3");
	}
#ifdef CUDA
	static int semid = -1;
	struct sembuf op_lock[2]=  { 0, 0, 0, /* wait for sem #0 to become 0 */
	     						 0, 1, SEM_UNDO /* then increment sem #0 by 1 */ };
	struct sembuf op_unlock[1]= { 0, -1, (IPC_NOWAIT | SEM_UNDO) /* decrement sem #0 by 1 (sets it to 0) */ };

	size_t mem_req, ws_sz;
	int Faux_sz = padoridim*(padoridim/2+1);
	int n[3] = {padoridim, padoridim, padoridim};
	cufftType cufft_type = CUFFT_R2C;
	if(sizeof(RFLOAT) == sizeof(double))
		cufft_type = CUFFT_D2Z;
	
	if(ref_dim == 3)
		Faux_sz *= padoridim;

	mem_req =  (size_t)1024;
	if(do_heavy && do_gpu)
	{
	    cufftEstimateMany(ref_dim, n, NULL, 0, 0, NULL, 0, 0, cufft_type, 1, &ws_sz);

		mem_req = (size_t)sizeof(RFLOAT)*MULTIDIM_SIZE(vol_in) +                   // dvol
				  (size_t)sizeof(Complex)*Faux_sz +                                // dFaux
				  (size_t)sizeof(RFLOAT)*MULTIDIM_SIZE(Mpad) +                     // dMpad
				  ws_sz + 4096;                                                    // workspace for cuFFT + extra space for alingment
	}

	CudaCustomAllocator *allocator = NULL;
	if(do_gpu && do_heavy)
	{
		int devid;
		size_t mem_free, mem_tot;
		cudaDeviceProp devProp;
		if (semid <0) 
		{
			HANDLE_ERROR(cudaGetDevice(&devid));
			if ( ( semid=semget(SEMKEY+devid, 1, IPC_CREAT | PERMS )) < 0 ) 
				REPORT_ERROR("semget error"); 
		}
		if (semop(semid, &op_lock[0], 2) < 0) 
			REPORT_ERROR("semop lock error");

		HANDLE_ERROR(cudaMemGetInfo(&mem_free, &mem_tot));
		if(mem_free > mem_req)
			allocator = new CudaCustomAllocator(mem_req, (size_t)16);
		else
		{
			do_gpu = false; // change local copy of do_gpu variable
			if (semop(semid, &op_unlock[0], 1) < 0)
				REPORT_ERROR("semop unlock error");
		}
	}
	AccPtrFactory ptrFactory(allocator);
	AccPtr<RFLOAT> dMpad = ptrFactory.make<RFLOAT>(MULTIDIM_SIZE(Mpad));
	AccPtr<Complex> dFaux = ptrFactory.make<Complex>(Faux_sz);
	AccPtr<RFLOAT> dvol = ptrFactory.make<RFLOAT>(MULTIDIM_SIZE(vol_in));
	if(do_heavy && do_gpu)
	{
		dvol.setHostPtr(MULTIDIM_ARRAY(vol_in));
		dvol.accAlloc();
		dvol.cpToDevice();
	}
#endif
	TIMING_TOC(TIMING_INIT1);

	TIMING_TIC(TIMING_GRID);
	// First do a gridding pre-correction on the real-space map:
	// Divide by the inverse Fourier transform of the interpolator in Fourier-space
	// 10feb11: at least in 2D case, this seems to be the wrong thing to do!!!
	// TODO: check what is best for subtomo!
	if (do_gridding)// && data_dim != 3)
	{
		if(do_heavy)
#ifdef CUDA
			if(do_gpu)
			{
				vol_in.setXmippOrigin();
				run_griddingCorrect(~dvol, interpolator, (RFLOAT)(ori_size * padding_factor), r_min_nn,
							    XSIZE(vol_in), YSIZE(vol_in), ZSIZE(vol_in));
			}
			else
#endif
			griddingCorrect(vol_in);
		else
			vol_in.setXmippOrigin();
	}

	TIMING_TOC(TIMING_GRID);

	TIMING_TIC(TIMING_PAD);

	// Pad translated map with zeros
	vol_in.setXmippOrigin();
	Mpad.setXmippOrigin();
	if(do_heavy)
	{
#ifdef CUDA
		if(do_gpu)
		{
			dMpad.accAlloc();
			run_padTranslatedMap(~dvol, ~dMpad,
				STARTINGX(vol_in),FINISHINGX(vol_in),STARTINGY(vol_in),FINISHINGY(vol_in),STARTINGZ(vol_in),FINISHINGZ(vol_in),   //Input dimensions
				STARTINGX(Mpad),  FINISHINGX(Mpad),  STARTINGY(Mpad),  FINISHINGY(Mpad),  STARTINGZ(Mpad),  FINISHINGZ(Mpad)      //Output dimensions
				);
			dvol.freeDevice();
		}
		else
#endif
		FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in) // This will also work for 2D
			A3D_ELEM(Mpad, k, i, j) = A3D_ELEM(vol_in, k, i, j);
	}
	TIMING_TOC(TIMING_PAD);

	TIMING_TIC(TIMING_TRANS);
	// Calculate the oversampled Fourier transform
	if(do_heavy)
#ifdef CUDA
		if(do_gpu)
		{
			dFaux.accAlloc();
			cufftResult err;
			cufftHandle plan;
			
			err = cufftCreate(&plan);
			if(err != CUFFT_SUCCESS)
				REPORT_ERROR("failed to create cufft plan");
			cufftSetAutoAllocation(plan, 0); // do not allocate work area
			//Allocate space with smart allocator
			AccPtr<char> fft_ws = ptrFactory.make<char>(ws_sz);
			fft_ws.accAlloc();
			cufftSetWorkArea(plan, ~fft_ws);
			err = cufftMakePlanMany(plan, ref_dim, n, NULL, 0, 0, NULL, 0, 0, cufft_type, 1, &ws_sz);
			if(err != CUFFT_SUCCESS)
				REPORT_ERROR("failed to create cufft plan");
			
			// do inverse FFT (dMpad->dFaux)
			if(sizeof(RFLOAT) == sizeof(double))
				err = cufftExecD2Z(plan,(cufftDoubleReal*)~dMpad, (cufftDoubleComplex*)~dFaux);
			else
				err = cufftExecR2C(plan,(cufftReal*)~dMpad, (cufftComplex*)~dFaux);
			if(err != CUFFT_SUCCESS)
				REPORT_ERROR("failed to exec fft");
			// deallocate plan, free mem
			cufftDestroy(plan);
			fft_ws.freeIfSet();
			
			size_t normfft = (size_t)padoridim*(size_t)padoridim;
			if(ref_dim == 3) normfft *= (size_t)padoridim;

			if(ref_dim == 2) Faux.reshape(padoridim,(padoridim/2+1));
			if(ref_dim == 3) Faux.reshape(padoridim,padoridim,(padoridim/2+1));
			
			scale((RFLOAT*)~dFaux, 2*dFaux.getSize(), 1./(RFLOAT)normfft);
		}
		else
#endif
		transformer.FourierTransform(Mpad, Faux, false);
	TIMING_TOC(TIMING_TRANS);

	TIMING_TIC(TIMING_CENTER);
	// Translate padded map to put origin of FT in the center
	if(do_heavy)
#ifdef CUDA
		if(do_gpu)
		{
			run_CenterFFTbySign(~dFaux, XSIZE(Faux), YSIZE(Faux), ZSIZE(Faux));
		}
		else
#endif
		CenterFFTbySign(Faux);
	TIMING_TOC(TIMING_CENTER);

	TIMING_TIC(TIMING_INIT2);
	// Free memory: Mpad no longer needed
#ifdef CUDA
	dMpad.freeIfSet();
#endif
	Mpad.clear();

	// Resize data array to the right size and initialise to zero
	initZeros(current_size);

	// Fill data only for those points with distance to origin less than max_r
	// (other points will be zero because of initZeros() call above
	// Also calculate radial power spectrum
#ifdef CUDA
	int fourier_mask_sz = (do_fourier_mask)?MULTIDIM_SIZE(*fourier_mask):16;
	int fmXsz, fmYsz, fmZsz;
	AccPtr<Complex> ddata = ptrFactory.make<Complex>(MULTIDIM_SIZE(data));
	AccPtr<RFLOAT> dfourier_mask = ptrFactory.make<RFLOAT>(fourier_mask_sz);
	AccPtr<RFLOAT> dpower_spectrum = ptrFactory.make<RFLOAT>(ori_size / 2 + 1);
	AccPtr<RFLOAT> dcounter = ptrFactory.make<RFLOAT>(ori_size / 2 + 1);

	if(do_heavy && do_gpu){
		ddata.accAlloc();
		dpower_spectrum.accAlloc();
		dcounter.accAlloc();
		ddata.deviceInit(0);
		dpower_spectrum.deviceInit(0);
		dcounter.deviceInit(0);
		dfourier_mask.accAlloc();
		fmXsz = fmYsz = fmZsz = 0;
		if(do_fourier_mask)
		{
			dfourier_mask.setHostPtr(MULTIDIM_ARRAY(*fourier_mask));
			dfourier_mask.cpToDevice();
			fmXsz = XSIZE(*fourier_mask);
			fmYsz = YSIZE(*fourier_mask);
			fmZsz = ZSIZE(*fourier_mask);
		}
	}
#endif
	power_spectrum.initZeros(ori_size / 2 + 1);
	MultidimArray<RFLOAT> counter(power_spectrum);
	counter.initZeros();
	TIMING_TOC(TIMING_INIT2);

	TIMING_TIC(TIMING_FAUX);
	int max_r2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
	int min_r2 = -1;
	if (min_ires > 0)
	{
		min_r2 = ROUND(min_ires * padding_factor) * ROUND(min_ires * padding_factor);
	}

	if(do_heavy)
	{
		RFLOAT weight = 1.;
#ifdef CUDA
		if(do_gpu)
		{
			run_calcPowerSpectrum(~dFaux, padoridim, ~ddata, YSIZE(data), ~dpower_spectrum, ~dcounter,
								  max_r2, min_r2, normfft, padding_factor, weight,
								  ~dfourier_mask, fmXsz, fmYsz, fmZsz, do_fourier_mask, ref_dim == 3);
			ddata.setHostPtr(MULTIDIM_ARRAY(data));
			ddata.cpToHost();
			dfourier_mask.freeIfSet();
		}
		else
#endif
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) // This will also work for 2D
		{
			int r2 = kp*kp + ip*ip + jp*jp;
			// The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
			// Set data array
			if (r2 <= max_r2)
			{
				if (do_fourier_mask) weight = FFTW_ELEM(*fourier_mask, ROUND(kp/padding_factor), ROUND(ip/padding_factor), ROUND(jp/padding_factor));
				// Set data array
				A3D_ELEM(data, kp, ip, jp) = weight * DIRECT_A3D_ELEM(Faux, k, i, j) * normfft;

				// Calculate power spectrum
				int ires = ROUND( sqrt((RFLOAT)r2) / padding_factor );
				// Factor two because of two-dimensionality of the complex plane
				DIRECT_A1D_ELEM(power_spectrum, ires) += norm(A3D_ELEM(data, kp, ip, jp)) / 2.;
				DIRECT_A1D_ELEM(counter, ires) += weight;

				// Apply high pass filter of the reference only after calculating the power spectrum
				if (r2 <= min_r2)
					A3D_ELEM(data, kp, ip, jp) = 0;
			}
		}
	}
	TIMING_TOC(TIMING_FAUX);

	/*
	FourierTransformer ft2;
	MultidimArray<Complex> Faux2(padding_factor * ori_size, (padding_factor * ori_size)/2+1);
	Image<RFLOAT> tt2(padding_factor * ori_size, padding_factor * ori_size);
	decenter(data, Faux2, max_r2);
	CenterFFTbySign(Faux2);
	windowFourierTransform(Faux2, padding_factor * ori_size);
	ft2.inverseFourierTransform(Faux2, tt2());
	tt2().setXmippOrigin();
	tt2().window(FIRST_XMIPP_INDEX(ori_size), FIRST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size), LAST_XMIPP_INDEX(ori_size));
	tt2.write("Fdata_proj.spi");
	std::cerr << "written Fdata_proj.spi" << std::endl;
	REPORT_ERROR("STOP");
	*/

	TIMING_TIC(TIMING_POW);
	// Calculate radial average of power spectrum
	if(do_heavy)
	{
#ifdef CUDA
		if(do_gpu)
		{
			run_updatePowerSpectrum(~dcounter, dcounter.getSize(), ~dpower_spectrum);
			dpower_spectrum.setHostPtr(MULTIDIM_ARRAY(power_spectrum));
			dpower_spectrum.cpToHost();
		}
		else
#endif
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(power_spectrum)
		{
			if (DIRECT_A1D_ELEM(counter, i) < 1.)
				DIRECT_A1D_ELEM(power_spectrum, i) = 0.;
			else
				DIRECT_A1D_ELEM(power_spectrum, i) /= DIRECT_A1D_ELEM(counter, i);
		}
	}
	TIMING_TOC(TIMING_POW);

	TIMING_TOC(TIMING_TOP);
#ifdef CUDA
	ddata.freeIfSet();
	dpower_spectrum.freeIfSet();
	dcounter.freeIfSet();
	dvol.freeIfSet();
	dMpad.freeIfSet();
	dFaux.freeIfSet();

	if(allocator != NULL)
	{
		delete allocator;
		if (semop(semid, &op_unlock[0], 1) < 0)
			REPORT_ERROR("semop unlock error");
	}
#endif

#ifdef PROJ_TIMING
	proj_timer.printTimes(false);
#endif

}

void Projector::griddingCorrect(MultidimArray<RFLOAT> &vol_in)
{
	// Correct real-space map by dividing it by the Fourier transform of the interpolator(s)
	vol_in.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in)
	{
		RFLOAT r = sqrt((RFLOAT)(k*k+i*i+j*j));
		// if r==0: do nothing (i.e. divide by 1)
		if (r > 0.)
		{
			RFLOAT rval = r / (ori_size * padding_factor);
			RFLOAT sinc = sin(PI * rval) / ( PI * rval);
			//RFLOAT ftblob = blob_Fourier_val(rval, blob) / blob_Fourier_val(0., blob);
			// Interpolation (goes with "interpolator") to go from arbitrary to fine grid
			if (interpolator==NEAREST_NEIGHBOUR && r_min_nn == 0)
			{
				// NN interpolation is convolution with a rectangular pulse, which FT is a sinc function
				A3D_ELEM(vol_in, k, i, j) /= sinc;
			}
			else if (interpolator==TRILINEAR || (interpolator==NEAREST_NEIGHBOUR && r_min_nn > 0) )
			{
				// trilinear interpolation is convolution with a triangular pulse, which FT is a sinc^2 function
				A3D_ELEM(vol_in, k, i, j) /= sinc * sinc;
			}
			else
				REPORT_ERROR("BUG Projector::griddingCorrect: unrecognised interpolator scheme.");
//#define DEBUG_GRIDDING_CORRECT
#ifdef DEBUG_GRIDDING_CORRECT
			if (k==0 && i==0 && j > 0)
				std::cerr << " j= " << j << " sinc= " << sinc << std::endl;
#endif
		}
	}
}

void Projector::project(MultidimArray<Complex > &f2d, Matrix2D<RFLOAT> &A)
{
	// f2d should already be in the right size (ori_size,orihalfdim)
	// AND the points outside r_max should already be zero...
	// f2d.initZeros();

	// Use the inverse matrix

	Matrix2D<RFLOAT> Ainv = A.inv();
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	// The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
	const int r_max_out = XSIZE(f2d) - 1;
	const int r_max_out_2 = r_max_out * r_max_out;

	const int r_max_ref = r_max * padding_factor;
	const int r_max_ref_2 = r_max_ref * r_max_ref;

	const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

//#define DEBUG
#ifdef DEBUG
	std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
	std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
	std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
	std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
	std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
	std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
	std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
	std::cerr << " max_r= "<< r_max << std::endl;
	std::cerr << " Ainv= " << Ainv << std::endl;
#endif

	for (int i = 0; i < YSIZE(f2d); i++)
	{
		const int y = (i <= r_max_out)? i : i - YSIZE(f2d);
		const int y2 = y * y;

		const int x_max = FLOOR(sqrt(r_max_out_2 - y2));

		for (int x = 0; x <= x_max; x++)
		{
			// sqrt(x*x + y*y) guaranteed to be < r_max_out

			// Get logical coordinates in the 3D map
			RFLOAT xp = Ainv(0,0) * x + Ainv(0,1) * y;
			RFLOAT yp = Ainv(1,0) * x + Ainv(1,1) * y;
			RFLOAT zp = Ainv(2,0) * x + Ainv(2,1) * y;

			const RFLOAT r_ref_2 = xp*xp + yp*yp + zp*zp;

			if (r_ref_2 > r_max_ref_2) continue;

			if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2)
			{
				// Only asymmetric half is stored

				bool is_neg_x = (xp < 0);

				if (is_neg_x)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					zp = -zp;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
				const int x0 = FLOOR(xp);
				const RFLOAT fx = xp - x0;
				const int x1 = x0 + 1;

				int y0 = FLOOR(yp);
				const RFLOAT fy = yp - y0;
				y0 -=  STARTINGY(data);
				const int y1 = y0 + 1;

				int z0 = FLOOR(zp);
				const RFLOAT fz = zp - z0;
				z0 -= STARTINGZ(data);
				const int z1 = z0 + 1;

				// Avoid reading outside the box
				if (x0 < 0 || x0+1 >= data.xdim
				 || y0 < 0 || y0+1 >= data.ydim
				 || z0 < 0 || z0+1 >= data.zdim)
				{
					continue;
				}

				// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
				const Complex d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
				const Complex d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
				const Complex d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
				const Complex d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
				const Complex d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
				const Complex d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
				const Complex d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
				const Complex d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

				// Set the interpolated value in the 2D output array
				const Complex dx00 = LIN_INTERP(fx, d000, d001);
				const Complex dx01 = LIN_INTERP(fx, d100, d101);
				const Complex dx10 = LIN_INTERP(fx, d010, d011);
				const Complex dx11 = LIN_INTERP(fx, d110, d111);
				const Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
				const Complex dxy1 = LIN_INTERP(fy, dx01, dx11);

				DIRECT_A2D_ELEM(f2d, i, x) = LIN_INTERP(fz, dxy0, dxy1);

				// Take complex conjugated for half with negative x
				if (is_neg_x)
				{
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));
				}

			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR ) // never actually used
			{
				int x0 = ROUND(xp);
				int y0 = ROUND(yp);
				int z0 = ROUND(zp);

				const bool is_neg_x = (x0 < 0);

				if (is_neg_x)
				{
					// Get complex conjugated hermitian symmetry pair
					x0 = -x0;
					y0 = -y0;
					z0 = -z0;
				}

				const int xr = x0 - STARTINGX(data);
				const int yr = y0 - STARTINGY(data);
				const int zr = z0 - STARTINGZ(data);

				if (xr < 0 || xr >= data.xdim
				 || yr < 0 || yr >= data.ydim
				 || zr < 0 || zr >= data.zdim)
				{
					continue;
				}

				if (is_neg_x)
				{
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A3D_ELEM(data, zr, yr, xr));
				}
				else
				{
					DIRECT_A2D_ELEM(f2d, i, x) = A3D_ELEM(data, zr, yr, xr);
				}

			} // endif NEAREST_NEIGHBOUR
			else
			{
				REPORT_ERROR("Unrecognized interpolator in Projector::project");
			}

		} // endif x-loop
	} // endif y-loop

#ifdef DEBUG
    std::cerr << "done with project..." << std::endl;
#endif
}

void Projector::projectGradient(Volume<t2Vector<Complex>>& img_out, Matrix2D<RFLOAT>& At)
{
	const int s = img_out.dimy;
	const int sh = img_out.dimx;

	Matrix2D<RFLOAT> Ainv = At.inv();

	// Go from the 2D slice coordinates to the 3D coordinates
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	for (int yy = 0; yy < s; yy++)
	{
		const double y = yy < sh? yy : yy - s;
		const double y2 = y * y;

		for (int xx = 0; xx < sh; xx++)
		{
			const double x = xx;

			if (x*x + y2 > sh*sh) continue;

			// Get logical coordinates in the 3D map
			double xp = Ainv(0,0) * x + Ainv(0,1) * y;
			double yp = Ainv(1,0) * x + Ainv(1,1) * y;
			double zp = Ainv(2,0) * x + Ainv(2,1) * y;

			bool is_neg_x;

			// Only asymmetric half is stored
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				is_neg_x = true;
			}
			else
			{
				is_neg_x = false;
			}

			// Trilinear interpolation (with physical coords)
			// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
			// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM

			int x0 = FLOOR(xp);
			double fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = FLOOR(yp);
			double fy = yp - y0;
			y0 -=  STARTINGY(data);
			int y1 = y0 + 1;

			int z0 = FLOOR(zp);
			double fz = zp - z0;
			z0 -= STARTINGZ(data);
			int z1 = z0 + 1;

			if (x0 < 0 || x0+1 >= data.xdim
			 || y0 < 0 || y0+1 >= data.ydim
			 || z0 < 0 || z0+1 >= data.zdim)
			{
				img_out(xx, yy, 0) = t2Vector<Complex>(Complex(0.0, 0.0), Complex(0.0, 0.0));
				continue;
			}

			Complex v000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
			Complex v001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
			Complex v010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
			Complex v011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
			Complex v100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
			Complex v101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
			Complex v110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
			Complex v111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

			Complex v00 = LIN_INTERP(fx, v000, v001);
			Complex v10 = LIN_INTERP(fx, v100, v101);
			Complex v01 = LIN_INTERP(fx, v010, v011);
			Complex v11 = LIN_INTERP(fx, v110, v111);

			Complex v0 = LIN_INTERP(fy, v00, v01);
			Complex v1 = LIN_INTERP(fy, v10, v11);

			// Complex v = LIN_INTERP(fz, v0, v1);

			Complex v00_dx = v001 - v000;
			Complex v10_dx = v101 - v100;
			Complex v01_dx = v011 - v010;
			Complex v11_dx = v111 - v110;
			Complex v0_dx = LIN_INTERP(fy, v00_dx, v01_dx);
			Complex v1_dx = LIN_INTERP(fy, v10_dx, v11_dx);
			Complex v_dx = LIN_INTERP(fz, v0_dx, v1_dx);

			Complex v0_dy = v01 - v00;
			Complex v1_dy = v11 - v10;
			Complex v_dy = LIN_INTERP(fz, v0_dy, v1_dy);

			Complex v_dz = v1 - v0;

			t3Vector<Complex> grad3D(v_dx, v_dy, v_dz);

			// Take complex conjugated for half with negative x
			if (is_neg_x)
			{
				grad3D.x = -(grad3D.x).conj();
				grad3D.y = -(grad3D.y).conj();
				grad3D.z = -(grad3D.z).conj();
			}

			img_out(xx, yy, 0).x = Ainv(0,0) * grad3D.x + Ainv(1,0) * grad3D.y + Ainv(2,0) * grad3D.z;
			img_out(xx, yy, 0).y = Ainv(0,1) * grad3D.x + Ainv(1,1) * grad3D.y + Ainv(2,1) * grad3D.z;
		} // endif x-loop
	} // endif y-loop
}

// Never actually used:
void Projector::project2Dto1D(MultidimArray<Complex > &f1d, Matrix2D<RFLOAT> &A)
{
	// f1d should already be in the right size (ori_size,orihalfdim)
	// AND the points outside r_max should already be zero...
	// f1d.initZeros();

	Matrix2D<RFLOAT> Ainv = A.inv();
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	// The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
	const int r_max_out = XSIZE(f1d) - 1;

	const int r_max_ref = r_max * padding_factor;
	const int r_max_ref_2 = r_max_ref * r_max_ref;

	const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;


	for (int x = 0; x <= r_max_out; x++)
	{
		// Get logical coordinates in the 2D map
		RFLOAT xp = Ainv(0,0) * x;
		RFLOAT yp = Ainv(1,0) * x;

		const RFLOAT r_ref_2 = xp*xp + yp*yp;

		if (r_ref_2 > r_max_ref_2) continue;

		if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2)
		{
			// Only asymmetric half is stored
			const bool is_neg_x = xp < 0;

			if (is_neg_x)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
			}

			// Trilinear interpolation (with physical coords)
			// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
			// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
			const int x0 = FLOOR(xp);
			const RFLOAT fx = xp - x0;
			const int x1 = x0 + 1;

			int y0 = FLOOR(yp);
			const RFLOAT fy = yp - y0;
			y0 -=  STARTINGY(data);
			const int y1 = y0 + 1;

			// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
			const Complex d00 = DIRECT_A2D_ELEM(data, y0, x0);
			const Complex d01 = DIRECT_A2D_ELEM(data, y0, x1);
			const Complex d10 = DIRECT_A2D_ELEM(data, y1, x0);
			const Complex d11 = DIRECT_A2D_ELEM(data, y1, x1);

			// Set the interpolated value in the 2D output array
			const Complex dx0 = LIN_INTERP(fx, d00, d01);
			const Complex dx1 = LIN_INTERP(fx, d10, d11);

			DIRECT_A1D_ELEM(f1d, x) = LIN_INTERP(fy, dx0, dx1);

			// Take complex conjugated for half with negative x
			if (is_neg_x)
			{
				DIRECT_A1D_ELEM(f1d, x) = conj(DIRECT_A1D_ELEM(f1d, x));
			}

		} // endif TRILINEAR
		else if (interpolator == NEAREST_NEIGHBOUR )
		{
			const int x0 = ROUND(xp);
			const int y0 = ROUND(yp);

			if (x0 < 0)
			{
				DIRECT_A1D_ELEM(f1d, x) = conj(A2D_ELEM(data, -y0, -x0));
			}
			else
			{
				DIRECT_A1D_ELEM(f1d, x) = A2D_ELEM(data, y0, x0);
			}

		} // endif NEAREST_NEIGHBOUR
		else
		{
			REPORT_ERROR("Unrecognized interpolator in Projector::project2Dto1D");
		}
	} // endif x-loop
}

void Projector::rotate2D(MultidimArray<Complex > &f2d, Matrix2D<RFLOAT> &A)
{
	// f2d should already be in the right size (ori_size,orihalfdim)
	// AND the points outside max_r should already be zero...
	// f2d.initZeros();

	// Use the inverse matrix
	Matrix2D<RFLOAT> Ainv = A.inv();
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	// The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
	const int r_max_out = XSIZE(f2d) - 1;
	const int r_max_out_2 = r_max_out * r_max_out;

	const int r_max_ref = r_max * padding_factor;
	const int r_max_ref_2 = r_max_ref * r_max_ref;

	const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

#ifdef DEBUG
	std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
	std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
	std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
	std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
	std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
	std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
	std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
	std::cerr << " max_r= "<< r_max << std::endl;
	std::cerr << " Ainv= " << Ainv << std::endl;
#endif

	for (int i=0; i < YSIZE(f2d); i++)
	{
		const int y = (i <= r_max_out)? i : i - YSIZE(f2d);
		const int y2 = y * y;

		const int x_max = FLOOR(sqrt(r_max_out_2 - y2));

		for (int x = 0; x <= x_max; x++)
		{
			// sqrt(x*x + y*y) guaranteed to be < r_max_out

			RFLOAT xp = Ainv(0,0) * x + Ainv(0,1) * y;
			RFLOAT yp = Ainv(1,0) * x + Ainv(1,1) * y;

			const int r_ref_2 = xp*xp + yp*yp;

			if (r_ref_2 > r_max_ref_2) continue;

			if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2)
			{
				const bool is_neg_x = xp < 0;

				// Only asymmetric half is stored
				if (is_neg_x)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
    			const int x0 = FLOOR(xp);
				const RFLOAT fx = xp - x0;
				const int x1 = x0 + 1;

				int y0 = FLOOR(yp);
				const RFLOAT fy = yp - y0;
				y0 -=  STARTINGY(data);
				const int y1 = y0 + 1;

				// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
				const Complex d00 = DIRECT_A2D_ELEM(data, y0, x0);
				const Complex d01 = DIRECT_A2D_ELEM(data, y0, x1);
				const Complex d10 = DIRECT_A2D_ELEM(data, y1, x0);
				const Complex d11 = DIRECT_A2D_ELEM(data, y1, x1);

				// Set the interpolated value in the 2D output array
				const Complex dx0 = LIN_INTERP(fx, d00, d01);
				const Complex dx1 = LIN_INTERP(fx, d10, d11);

				DIRECT_A2D_ELEM(f2d, i, x) = LIN_INTERP(fy, dx0, dx1);

				// Take complex conjugated for half with negative x
				if (is_neg_x)
				{
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));
				}
			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR ) // never used
			{
				const int x0 = ROUND(xp);
				const int y0 = ROUND(yp);

				if (x0 < 0)
				{
					DIRECT_A2D_ELEM(f2d, i, x) = conj(A2D_ELEM(data, -y0, -x0));
				}
				else
				{
					DIRECT_A2D_ELEM(f2d, i, x) = A2D_ELEM(data, y0, x0);
				}
			} // endif NEAREST_NEIGHBOUR
			else
			{
				REPORT_ERROR("Unrecognized interpolator in Projector::project");
			}
		} // endif x-loop
	} // endif y-loop
}


void Projector::rotate3D(MultidimArray<Complex > &f3d, Matrix2D<RFLOAT> &A)
{
	// f3d should already be in the right size (ori_size,orihalfdim)
	// AND the points outside max_r should already be zero
	// f3d.initZeros();

	// Use the inverse matrix
	Matrix2D<RFLOAT> Ainv = A.inv();
	Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

	const int r_max_out = XSIZE(f3d) - 1;
	const int r_max_out_2 = r_max_out * r_max_out;

	const int r_max_ref = r_max * padding_factor;
	const int r_max_ref_2 = r_max_ref * r_max_ref;

	const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

#ifdef DEBUG
	std::cerr << " XSIZE(f3d)= "<< XSIZE(f3d) << std::endl;
	std::cerr << " YSIZE(f3d)= "<< YSIZE(f3d) << std::endl;
	std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
	std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
	std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
	std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
	std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
	std::cerr << " max_r= "<< r_max << std::endl;
	std::cerr << " Ainv= " << Ainv << std::endl;
#endif

	for (int k = 0; k < ZSIZE(f3d); k++)
	{
		const int z = (k <= r_max_out)? k : k - ZSIZE(f3d);
		const int z2 = z * z;

		for (int i = 0; i < YSIZE(f3d); i++)
		{
			const int y = (i <= r_max_out)? i : i - YSIZE(f3d);
			const int y2 = y * y;

			const RFLOAT yz2 = y2 + z2;

			// avoid negative square root
			if (yz2 > r_max_out_2) continue;

			const int x_max = FLOOR(sqrt(r_max_out_2 - yz2));

			for (int x = 0; x <= x_max; x++)
			{
				// Get logical coordinates in the 3D map
				RFLOAT xp = Ainv(0,0) * x + Ainv(0,1) * y + Ainv(0,2) * z;
				RFLOAT yp = Ainv(1,0) * x + Ainv(1,1) * y + Ainv(1,2) * z;
				RFLOAT zp = Ainv(2,0) * x + Ainv(2,1) * y + Ainv(2,2) * z;

				const int r_ref_2 = xp*xp + yp*yp + zp*zp;

				if (r_ref_2 > r_max_ref_2) continue;

				if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2)
				{
					// Only asymmetric half is stored
					const bool is_neg_x = (xp < 0);

					if (is_neg_x)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
					}

					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					const int x0 = FLOOR(xp);
					const RFLOAT fx = xp - x0;
					const int x1 = x0 + 1;

					int y0 = FLOOR(yp);
					const RFLOAT fy = yp - y0;
					y0 -=  STARTINGY(data);
					const int y1 = y0 + 1;

					int z0 = FLOOR(zp);
					const RFLOAT fz = zp - z0;
					z0 -=  STARTINGZ(data);
					const int z1 = z0 + 1;

					// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
					const Complex d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
					const Complex d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
					const Complex d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
					const Complex d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
					const Complex d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
					const Complex d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
					const Complex d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
					const Complex d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

					// Set the interpolated value in the 2D output array
					// interpolate in x
					const Complex dx00 = LIN_INTERP(fx, d000, d001);
					const Complex dx01 = LIN_INTERP(fx, d100, d101);
					const Complex dx10 = LIN_INTERP(fx, d010, d011);
					const Complex dx11 = LIN_INTERP(fx, d110, d111);
					// interpolate in y
					const Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
					const Complex dxy1 = LIN_INTERP(fy, dx01, dx11);
					//interpolate in z
					DIRECT_A3D_ELEM(f3d, k, i, x) = LIN_INTERP(fz, dxy0, dxy1);

					// Take complex conjugated for half with negative x
					if (is_neg_x)
					{
						DIRECT_A3D_ELEM(f3d, k, i, x) = conj(DIRECT_A3D_ELEM(f3d, k, i, x));
					}

				} // endif TRILINEAR
				else if (interpolator == NEAREST_NEIGHBOUR )
				{
					const int x0 = ROUND(xp);
					const int y0 = ROUND(yp);
					const int z0 = ROUND(zp);

					if (x0 < 0)
					{
						DIRECT_A3D_ELEM(f3d, k, i, x) = conj(A3D_ELEM(data, -z0, -y0, -x0));
					}
					else
					{
						DIRECT_A3D_ELEM(f3d, k, i, x) = A3D_ELEM(data, z0, y0, x0);
					}

				} // endif NEAREST_NEIGHBOUR
				else
				{
					REPORT_ERROR("Unrecognized interpolator in Projector::project");
				}
			} // endif x-loop
		} // endif y-loop
	} // endif z-loop
}

