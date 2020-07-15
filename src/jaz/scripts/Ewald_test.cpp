#include <src/args.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/projection/point_insertion.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/optics/ctf_helper.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/math/Euler_angles.h>
#include <src/CPlot2D.h>
#include <src/backprojector.h>

using namespace gravis;

#define INPUT_PRECISION float
#define OUTPUT_PRECISION float

int main(int argc, char *argv[])
{
	
	const int s = 180;
	const int sh = s/2 + 1;
	const int num_observations = 1;		
	const double SNR = 1.0;
	const double Ewald_radius = 2 * s;
	
	bool forward, slow, explicit_gridding, legacy_backprojector, crossed;	
	int num_threads;
	
	
	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		forward = parser.checkOption("--forward", "Use forward mapping");
		slow = parser.checkOption("--slow", "Use slow backward mapping");
		explicit_gridding = parser.checkOption("--xg", "Use explicit gridding");
		legacy_backprojector = parser.checkOption("--old", "Use legacy backprojector");
		
		num_threads = textToInteger(parser.getOption("--j", "Number of threads", "6"));
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	
	const std::string tag 
		= "Ewald_" + (
	            legacy_backprojector? 
				"legacy" 
	              :
				(   
					std::string("new_")
				  + (forward? "forward" : "backward") 
				  + (slow? "_slow" : "")
				  + (crossed? "_crossed" : "")
				));
	
	std::cout << "tag: " << tag << std::endl;
		
	BufferedImage<INPUT_PRECISION> observation_P_RS(s,s);
	BufferedImage<INPUT_PRECISION> observation_Q_RS(s,s);
	
	const double falloff_sigma2 = s*s/64.0;
	        
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double xx = x < s/2? x : x - s;
		const double yy = y < s/2? y : y - s;
		const double r2 = xx*xx + yy*yy;
		
		observation_P_RS(x,y) = 2 * exp(-0.5*r2/falloff_sigma2) * (rand() / (double)RAND_MAX + 0.5);
		observation_Q_RS(x,y) = 2 * exp(-0.5*r2/falloff_sigma2) * (rand() / (double)RAND_MAX + 0.5);
	}
	
	observation_P_RS.write("observation_P_RS_"+tag+".mrc");
	observation_Q_RS.write("observation_Q_RS_"+tag+".mrc");
	
	BufferedImage<tComplex<INPUT_PRECISION>> observation_P, observation_Q;
	FFT::FourierTransform(observation_P_RS, observation_P, FFT::Both);
	FFT::FourierTransform(observation_Q_RS, observation_Q, FFT::Both);

	
	observation_P.writeVtk("observation_P_"+tag+".vtk");
	observation_Q.writeVtk("observation_Q_"+tag+".vtk");
	
	
	BufferedImage<INPUT_PRECISION> ctf(sh,s);
	ctf.fill(1);
	
	BufferedImage<tComplex<OUTPUT_PRECISION>> data(sh,s,s);
	BufferedImage<OUTPUT_PRECISION> weight(sh,s,s), spreading_function(sh,s,s);
	
	data.fill(tComplex<OUTPUT_PRECISION>(0,0));
	weight.fill(0);
	spreading_function.fill(0);
		
	BackProjector backprojector;
	Image<Complex> obervation_P_legacy, obervation_Q_legacy;
	Image<RFLOAT> ctf_legacy;
	
	if (legacy_backprojector)
	{
		observation_P.copyTo(obervation_P_legacy);
		observation_Q.copyTo(obervation_Q_legacy);
		ctf.copyTo(ctf_legacy);
		
		backprojector = BackProjector(s, 3, "C1", TRILINEAR, 1, 10, 0, 1.9, 15, 2, true);
		backprojector.initZeros(-1);
	}
	
	Log::beginProgress("Backprojecting observations", num_observations);
	
	for (int i = 0; i < num_observations; i++)
	{
		Log::updateProgress(i);
		
		double phi = 2.0 * PI * rand() / (double)RAND_MAX;
		double psi = 2.0 * PI * rand() / (double)RAND_MAX;		
		double tilt = (PI/2.0) * sin(PI * rand() / (double)RAND_MAX - PI/2.0);
		
		if (i == 0)
		{
			phi  = DEG2RAD(60);
			psi  = DEG2RAD(0);
			tilt = DEG2RAD(0);
		}
		
		d4Matrix proj = Euler::anglesToMatrix4(phi, tilt, psi);
				
		if (legacy_backprojector)
		{
			Matrix2D<RFLOAT> A(3,3);
			
			for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++)
			{
				A(r,c) = proj(r,c);
			}
			
			backprojector.set2DFourierTransform(
					obervation_P_legacy(), A, &ctf_legacy(), Ewald_radius, true);
			
			backprojector.set2DFourierTransform(
					obervation_Q_legacy(), A, &ctf_legacy(), Ewald_radius, false);
		}
		else
		{
			if (forward)
			{
				ClippedPointInsertion<INPUT_PRECISION,OUTPUT_PRECISION> clippedInsertion;
				
				FourierBackprojection::backprojectSphere_forward(
					clippedInsertion,
					observation_P, ctf, proj, Ewald_radius,
					data, weight);
				
				FourierBackprojection::backprojectSphere_forward(
					clippedInsertion,
					observation_Q, ctf, proj, -Ewald_radius,
					data, weight);
			}
			else
			{
				CtfHelper::CTFP_CTFQ_Pair<INPUT_PRECISION> ctfPQ 
				        = CtfHelper::stitchHalves(observation_P, observation_Q);
				
				ctfPQ.pq.writeVtk("observation_PQ_"+tag+".vtk");
				
				if (!slow)
				{
					FourierBackprojection::backprojectSphere_backward(
						ctfPQ.pq, ctf, proj, Ewald_radius,
						data, weight,
						num_threads);
					
					FourierBackprojection::backprojectSphere_backward(
						ctfPQ.qp, ctf, proj, -Ewald_radius,
						data, weight,
						num_threads);
				}
				else
				{
						
					FourierBackprojection::backprojectSphere_backward_slow(
						ctfPQ.pq, ctf, proj, Ewald_radius,
						data, weight,
						num_threads);
					
					FourierBackprojection::backprojectSphere_backward_slow(
						ctfPQ.qp, ctf, proj, -Ewald_radius,
						data, weight,
						num_threads);
				}
				
				
				if (explicit_gridding)
				{
					FourierBackprojection::backprojectSpreadingFunction(
						proj, spreading_function);
				}
			}
		}
	}
	
	Log::endProgress();
	
	
	BufferedImage<OUTPUT_PRECISION> data_div_RS(s,s,s);
	
	
	
	if (legacy_backprojector)
	{		
		MultidimArray<Complex> data_centered_xmipp(s,s,sh);
		backprojector.decenter(backprojector.data, data_centered_xmipp, s*s/4);
		RawImage<Complex> data_centered(data_centered_xmipp);
		data_centered.writeVtk("data_"+tag+".vtk");
		        
		Image<RFLOAT> vol_xmipp;		
		vol_xmipp().initZeros(s, s, s);
		vol_xmipp().setXmippOrigin();
		
		MultidimArray<RFLOAT> tau2;
		
		backprojector.reconstruct(vol_xmipp(), 10, false, tau2);
		
		data_div_RS.copyDataAndSizeFrom(vol_xmipp);
	}
	else
	{
		data.writeVtk("dataDebug_"+tag+".vtk");
		        
		BufferedImage<OUTPUT_PRECISION> data_RS(s,s,s);
		
		if (forward || !explicit_gridding)
		{
			Reconstruction::griddingCorrect3D_sinc2(
				data, data_RS, 
				true, 1);
		}
		else
		{
			Reconstruction::griddingCorrect3D(
				data, spreading_function,  // in
				data_RS,                   // out
				true, num_threads);
		}
		
		Reconstruction::ctfCorrect3D_Wiener(
			data_RS, weight,           // in
			data_div_RS,               // out
			1.0 / SNR, num_threads);
	}
	
	data_div_RS.write("reconstruction_"+tag+".mrc");
	
	BufferedImage<tComplex<OUTPUT_PRECISION>> data_div_FS;
	FFT::FourierTransform(data_div_RS, data_div_FS);
	data_div_FS.writeVtk("reconstruction_"+tag+"_FS.vtk");
	
	return 0;
}
