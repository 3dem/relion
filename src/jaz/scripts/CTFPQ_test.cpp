#include <src/image.h>
#include <src/ctf.h>
#include <src/mask.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/noise.h>
#include <src/jaz/util/zio.h>




void applyCTFPandCTFQ_original(
		MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
		MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, bool skip_mask,
        int newbox, int nr_sectors, bool is_reverse, 
        double angpix, double mask_diameter, double width_mask_edge)
{
	outP.resize(Fin);
	outQ.resize(Fin);
	float angle_step = 180./nr_sectors;
	
	for (float angle = 0.; angle < 180.;  angle +=angle_step)
	{
		MultidimArray<Complex> CTFP(Fin), Fapp(Fin);
		MultidimArray<RFLOAT> Iapp(YSIZE(Fin), YSIZE(Fin));
		
		// Two passes: one for CTFP, one for CTFQ
		for (int ipass = 0; ipass < 2; ipass++)
		{
			bool is_my_positive = (ipass == 1) ? is_reverse : !is_reverse;

			// Get CTFP and multiply the Fapp with it
			ctf.getCTFPImage(CTFP, YSIZE(Fin), YSIZE(Fin), angpix, is_my_positive, angle);

			Fapp = Fin * CTFP; // element-wise complex multiplication!

			if (!skip_mask)
			{
				// inverse transform and mask out the particle....
				CenterFFTbySign(Fapp);
				transformer.inverseFourierTransform(Fapp, Iapp);

				softMaskOutsideMap(Iapp, ROUND(mask_diameter/(angpix*2.)), (RFLOAT)width_mask_edge);

				// Re-box to a smaller size if necessary....
				if (newbox > 0 && newbox < YSIZE(Fin))
				{
					Iapp.setXmippOrigin();
					Iapp.window(FIRST_XMIPP_INDEX(newbox), FIRST_XMIPP_INDEX(newbox),
					            LAST_XMIPP_INDEX(newbox),  LAST_XMIPP_INDEX(newbox));
				}
				
				// Back into Fourier-space
				transformer.FourierTransform(Iapp, Fapp, false); // false means: leave Fapp in the transformer
				CenterFFTbySign(Fapp);
			}

			// First time round: resize the output arrays
			if (ipass == 0 && fabs(angle) < XMIPP_EQUAL_ACCURACY)
			{
				outP.resize(Fapp);
				outQ.resize(Fapp);
			}

			// Now set back the right parts into outP (first pass) or outQ (second pass)
			float anglemin = angle + 90. - (0.5*angle_step);
			float anglemax = angle + 90. + (0.5*angle_step);

			// angles larger than 180
			bool is_sector_reverse = false;
			
			if (anglemin >= 180.)
			{
				anglemin -= 180.;
				anglemax -= 180.;
				is_sector_reverse = true;
			}
			
			MultidimArray<Complex> *myCTFPorQ, *myCTFPorQb;
			
			if (is_sector_reverse)
			{
				myCTFPorQ  = (ipass == 0) ? &outQ : &outP;
				myCTFPorQb = (ipass == 0) ? &outP : &outQ;
			}
			else
			{
				myCTFPorQ  = (ipass == 0) ? &outP : &outQ;
				myCTFPorQb = (ipass == 0) ? &outQ : &outP;
			}

			// Deal with sectors with the Y-axis in the middle of the sector...
			bool do_wrap_max = false;
			if (anglemin < 180. && anglemax > 180.)
			{
				anglemax -= 180.;
				do_wrap_max = true;
			}

			// use radians instead of degrees
			anglemin = DEG2RAD(anglemin);
			anglemax = DEG2RAD(anglemax);
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(CTFP)
			{
				RFLOAT x = (RFLOAT)jp;
				RFLOAT y = (RFLOAT)ip;
				RFLOAT myangle = (x*x+y*y > 0) ? acos(y/sqrt(x*x+y*y)) : 0; // dot-product with Y-axis: (0,1)
				
				// Only take the relevant sector now...
				if (do_wrap_max)
				{
					if (myangle >= anglemin)
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
					else if (myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQb, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
				else
				{
					if (myangle >= anglemin && myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
			}
		}
	}
}



void applyCTFPandCTFQ_logged(
		MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
		MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, bool skip_mask,
        int newbox, int nr_sectors, bool is_reverse, 
        double angpix, double mask_diameter, double width_mask_edge)
{
	outP.resize(Fin);
	outQ.resize(Fin);
	float angle_step = 180./nr_sectors;
	
	outP.initZeros();
	outQ.initZeros();
	
	RawImage<Complex>(outP).writeVtk("initial_outP.vtk");
	RawImage<Complex>(outQ).writeVtk("initial_outQ.vtk");
	
	for (float angle = 0.; angle < 180.;  angle +=angle_step)
	{
		std::cout << "angle = " << angle << std::endl;
		
		
		MultidimArray<Complex> CTFP(Fin), Fapp(Fin);
		MultidimArray<RFLOAT> Iapp(YSIZE(Fin), YSIZE(Fin));
		
		// Two passes: one for CTFP, one for CTFQ
		for (int ipass = 0; ipass < 2; ipass++)
		{
			std::cout << "ipass = " << ipass << std::endl;
			
			bool is_my_positive = (ipass == 1) ? is_reverse : !is_reverse;

			// Get CTFP and multiply the Fapp with it
			ctf.getCTFPImage(CTFP, YSIZE(Fin), YSIZE(Fin), angpix, is_my_positive, angle);

			Fapp = Fin * CTFP; // element-wise complex multiplication!
			
			RawImage<Complex>(CTFP).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_CTFP.vtk");
			RawImage<Complex>(Fapp).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_Fapp_0.vtk");
			
			if (!skip_mask)
			{
				// inverse transform and mask out the particle....
				CenterFFTbySign(Fapp);
				transformer.inverseFourierTransform(Fapp, Iapp);

				softMaskOutsideMap(Iapp, ROUND(mask_diameter/(angpix*2.)), (RFLOAT)width_mask_edge);

				// Re-box to a smaller size if necessary....
				if (newbox > 0 && newbox < YSIZE(Fin))
				{
					Iapp.setXmippOrigin();
					Iapp.window(FIRST_XMIPP_INDEX(newbox), FIRST_XMIPP_INDEX(newbox),
					            LAST_XMIPP_INDEX(newbox),  LAST_XMIPP_INDEX(newbox));
				}
				
				// Back into Fourier-space
				transformer.FourierTransform(Iapp, Fapp, false); // false means: leave Fapp in the transformer
				CenterFFTbySign(Fapp);
			}
			
			RawImage<Complex>(Fapp).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_Fapp_1.vtk");

			// Now set back the right parts into outP (first pass) or outQ (second pass)
			float anglemin = angle + 90. - (0.5*angle_step);
			float anglemax = angle + 90. + (0.5*angle_step);

			// angles larger than 180
			bool is_sector_reverse = false;
			
			if (anglemin >= 180.)
			{
				anglemin -= 180.;
				anglemax -= 180.;
				is_sector_reverse = true;
			}
			
			std::cout << "\tanglemin = " << anglemin << std::endl;
			std::cout << "\tanglemax = " << anglemax << std::endl;
			std::cout << "\tis_sector_reverse = " << is_sector_reverse << std::endl;
			
			MultidimArray<Complex> *myCTFPorQ, *myCTFPorQb;
			
			if (is_sector_reverse)
			{
				myCTFPorQ  = (ipass == 0) ? &outQ : &outP;
				myCTFPorQb = (ipass == 0) ? &outP : &outQ;
			}
			else
			{
				myCTFPorQ  = (ipass == 0) ? &outP : &outQ;
				myCTFPorQb = (ipass == 0) ? &outQ : &outP;
			}

			// Deal with sectors with the Y-axis in the middle of the sector...
			bool do_wrap_max = false;
			if (anglemin < 180. && anglemax > 180.)
			{
				anglemax -= 180.;
				do_wrap_max = true;
			}
			
			std::cout << "\tdo_wrap_max = " << do_wrap_max << std::endl;

			// use radians instead of degrees
			anglemin = DEG2RAD(anglemin);
			anglemax = DEG2RAD(anglemax);
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(CTFP)
			{
				RFLOAT x = (RFLOAT)jp;
				RFLOAT y = (RFLOAT)ip;
				RFLOAT myangle = (x*x+y*y > 0) ? acos(y/sqrt(x*x+y*y)) : 0; // dot-product with Y-axis: (0,1)
				
				// Only take the relevant sector now...
				if (do_wrap_max)
				{
					if (myangle >= anglemin)
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
					else if (myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQb, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
				else
				{
					if (myangle >= anglemin && myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
			}
			
			RawImage<Complex>(outP).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_outP.vtk");
			RawImage<Complex>(outQ).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_outQ.vtk");
		}
	}
}

void applyCTFPandCTFQ_refac(
		MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
		MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, bool skip_mask,
        int newbox, int nr_sectors, bool is_reverse, 
        double angpix, double mask_diameter, double width_mask_edge)
{
	outP.resize(Fin);
	outQ.resize(Fin);
	float angle_step = 180./nr_sectors;
	
	outP.initZeros();
	outQ.initZeros();
	
	RawImage<Complex>(outP).writeVtk("initial_outP.vtk");
	RawImage<Complex>(outQ).writeVtk("initial_outQ.vtk");
	
	for (float angle = 0.; angle < 180.;  angle +=angle_step)
	{
		std::cout << "angle = " << angle << std::endl;
		
		
		MultidimArray<Complex> CTFP(Fin), Fapp(Fin);
		MultidimArray<RFLOAT> Iapp(YSIZE(Fin), YSIZE(Fin));
		
		// Two passes: one for CTFP, one for CTFQ
		for (int ipass = 0; ipass < 2; ipass++)
		{
			std::cout << "ipass = " << ipass << std::endl;
			
			bool is_my_positive = (ipass == 1) ? is_reverse : !is_reverse;

			// Get CTFP and multiply the Fapp with it
			ctf.getCTFPImage(CTFP, YSIZE(Fin), YSIZE(Fin), angpix, is_my_positive, angle);

			Fapp = Fin * CTFP; // element-wise complex multiplication!
			
			RawImage<Complex>(CTFP).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_CTFP.vtk");
			RawImage<Complex>(Fapp).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_Fapp_0.vtk");
			
			if (!skip_mask)
			{
				// inverse transform and mask out the particle....
				CenterFFTbySign(Fapp);
				transformer.inverseFourierTransform(Fapp, Iapp);

				softMaskOutsideMap(Iapp, ROUND(mask_diameter/(angpix*2.)), (RFLOAT)width_mask_edge);

				// Re-box to a smaller size if necessary....
				if (newbox > 0 && newbox < YSIZE(Fin))
				{
					Iapp.setXmippOrigin();
					Iapp.window(FIRST_XMIPP_INDEX(newbox), FIRST_XMIPP_INDEX(newbox),
					            LAST_XMIPP_INDEX(newbox),  LAST_XMIPP_INDEX(newbox));
				}
				
				// Back into Fourier-space
				transformer.FourierTransform(Iapp, Fapp, false); // false means: leave Fapp in the transformer
				CenterFFTbySign(Fapp);
			}
			
			RawImage<Complex>(Fapp).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_Fapp_1.vtk");

			// Now set back the right parts into outP (first pass) or outQ (second pass)
			float anglemin = angle + 90. - (0.5*angle_step);
			float anglemax = angle + 90. + (0.5*angle_step);

			// angles larger than 180
			bool is_sector_reverse = false;
			
			if (anglemin >= 180.)
			{
				anglemin -= 180.;
				anglemax -= 180.;
				is_sector_reverse = true;
			}
			
			std::cout << "\tanglemin = " << anglemin << std::endl;
			std::cout << "\tanglemax = " << anglemax << std::endl;
			std::cout << "\tis_sector_reverse = " << is_sector_reverse << std::endl;
			
			MultidimArray<Complex> *myCTFPorQ, *myCTFPorQb;
			
			if (is_sector_reverse)
			{
				myCTFPorQ  = (ipass == 0) ? &outQ : &outP;
				myCTFPorQb = (ipass == 0) ? &outP : &outQ;
			}
			else
			{
				myCTFPorQ  = (ipass == 0) ? &outP : &outQ;
				myCTFPorQb = (ipass == 0) ? &outQ : &outP;
			}

			// Deal with sectors with the Y-axis in the middle of the sector...
			bool do_wrap_max = false;
			if (anglemin < 180. && anglemax > 180.)
			{
				anglemax -= 180.;
				do_wrap_max = true;
			}
			
			std::cout << "\tdo_wrap_max = " << do_wrap_max << std::endl;

			// use radians instead of degrees
			anglemin = DEG2RAD(anglemin);
			anglemax = DEG2RAD(anglemax);
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(CTFP)
			{
				RFLOAT x = (RFLOAT)jp;
				RFLOAT y = (RFLOAT)ip;
				RFLOAT myangle = (x*x+y*y > 0) ? acos(y/sqrt(x*x+y*y)) : 0; // dot-product with Y-axis: (0,1)
				
				// Only take the relevant sector now...
				if (do_wrap_max)
				{
					if (myangle >= anglemin)
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
					else if (myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQb, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
				else
				{
					if (myangle >= anglemin && myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
			}
			
			RawImage<Complex>(outP).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_outP.vtk");
			RawImage<Complex>(outQ).writeVtk(ZIO::itoa(angle)+"_"+ZIO::itoa(ipass)+"_outQ.vtk");
		}
	}
}


int main()
{
	const int s = 180;
	const int sh = s/2 + 1;
	
	const bool skip_mask = false;
	const int newbox = s;
	const int nr_sectors = 2;
	const bool is_reverse = false;
	const double angpix = 1;
	const double mask_diameter = s/2.0;
	const double width_mask_edge = 3.0;
		
	CTF ctf;
	ctf.setValues(8000, 12000, 30, 300, 2.7, 0.1, 0, 1, 0);
	
	
	MultidimArray<Complex> obs_xmipp(s,s,sh), obsP_xmipp, obsQ_xmipp;
	
	FourierTransformer transformer;
	
	BufferedImage<Complex> obs = ImageNoise::generateSmoothFourierNoise<RFLOAT>(s, s/4.0);
	
	obs.writeVtk("obs.vtk");
	obs.copyTo(obs_xmipp);
	
	applyCTFPandCTFQ_logged(
		obs_xmipp, ctf, transformer, obsP_xmipp, obsQ_xmipp, skip_mask,
		newbox, nr_sectors, is_reverse, angpix, mask_diameter, width_mask_edge);
	
	RawImage<Complex> obsP(obsP_xmipp);
	RawImage<Complex> obsQ(obsQ_xmipp);
	
	obs.writeVtk("obs.vtk");
	obsP.writeVtk("obsP.vtk");
	obsQ.writeVtk("obsQ.vtk");
	
	(obsP + obsQ).writeVtk("obsPplusQ.vtk");
	
	BufferedImage<Complex> flatCTF = obs;
	
	RFLOAT bs = (RFLOAT)s * angpix;
	
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const double xx = x / bs;		
		const double yy = (y < s/2? y : y - s) / bs;		
		
		flatCTF(x,y) *= ctf.getCTF(xx,yy);
	}
	
	flatCTF.writeVtk("obsCTF.vtk");
	
	
	MultidimArray<RFLOAT> Fctf(s,sh);
		
	ctf.applyWeightEwaldSphereCurvature_noAniso(
			Fctf, s, s, angpix, mask_diameter);
	
	RawImage<RFLOAT>(Fctf).write("maskedCTF.mrc");
	
	
	
	
	return 0;
}
