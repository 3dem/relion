#include <src/args.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
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
	const int num_shells = 5;
	const int num_spheres = 2 * num_shells - 1;
	const int num_threads = 6;
	const double outer_radius = s / 2;
	
	const std::string tag = "Ewald_"+ZIO::itoa(num_observations)+"_new_slow";
	
	const double SNR = 1.0;
	const double Ewald_radius = 5 * s;
	const bool forward = false;
	const bool slow = true;
	const bool explicit_gridding = false;
	const bool legacy_backprojector = false;
	const bool curved = true;
	
	std::vector<double> sphere_radius(num_spheres);
	std::vector<double> sphere_scale(num_spheres);
	
	for (int i = 0; i < num_spheres; i++)
	{
		const double r = (i + 1.0) * outer_radius / num_spheres;
		sphere_radius[i] = r;
		sphere_scale[i] = 1 - 2 * (i % 2);
	}
	
	BufferedImage<INPUT_PRECISION> observation(s,s);
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		double sum = 0.0;
		
		const double xx = x < s/2? x : x - s;
		const double yy = y < s/2? y : y - s;
		const double r2 = xx*xx + yy*yy;
		
		for (int i = 0; i < num_spheres; i++)
		{
			const double r2s = sphere_radius[i] * sphere_radius[i];
			const double d2 = r2s - r2;
			
			if (d2 > 0.0)
			{
				const double l = 2.0 * sqrt(d2);
				sum += sphere_scale[i] * l; 
			}
		}
		
		observation(x,y) = sum;
	}
	
	observation.write("observation_"+tag+".mrc");
	
	BufferedImage<tComplex<INPUT_PRECISION>> observation_FS;
	FFT::FourierTransform(observation, observation_FS, FFT::Both);
		
	BufferedImage<INPUT_PRECISION> ctf(sh,s);
	ctf.fill(1);
	
	BufferedImage<tComplex<OUTPUT_PRECISION>> data(sh,s,s);
	BufferedImage<OUTPUT_PRECISION> weight(sh,s,s), spreading_function(sh,s,s);
	
	data.fill(tComplex<OUTPUT_PRECISION>(0,0));
	weight.fill(0);
	spreading_function.fill(0);
		
	BackProjector backprojector;
	Image<Complex> obervation_FS_legacy;
	Image<RFLOAT> ctf_legacy;
	
	if (legacy_backprojector)
	{
		observation_FS.copyTo(obervation_FS_legacy);
		ctf.copyTo(ctf_legacy);
		
		backprojector = BackProjector(s, 3, "C1", TRILINEAR, 1, 10, 0, 1.9, 15, 2, true);
		backprojector.initZeros(-1);
	}
	
	Log::beginProgress("Backprojecting observations", num_observations);
	
	for (int i = 0; i < num_observations; i++)
	{
		Log::updateProgress(i);
		
		const double phi = 2.0 * PI * rand() / (double)RAND_MAX;
		const double psi = 2.0 * PI * rand() / (double)RAND_MAX;		
		const double tilt = (PI/2.0) * sin(PI * rand() / (double)RAND_MAX - PI/2.0);
		
		d4Matrix proj = Euler::anglesToMatrix4(phi, tilt, psi);
				
		if (legacy_backprojector)
		{
			Matrix2D<RFLOAT> A(3,3);
			
			for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++)
			{
				A(r,c) = proj(r,c);
			}
			
			backprojector.set2DFourierTransform(obervation_FS_legacy(), A, &ctf_legacy(), Ewald_radius);
		}
		else
		{
			if (forward)
			{
				/*if (wrap_voxels)
				{
					FourierBackprojection::backprojectSlice_noSF_fwd_wrap(
						observation_FS, ctf, proj,
						data, weight);
				}
				else
				{
					FourierBackprojection::backprojectSlice_noSF_fwd_clip(
						observation_FS, ctf, proj,
						data, weight);
				}*/
			}
			else
			{
				if (curved)
				{
					if (!slow)
					{
						FourierBackprojection::backprojectSlice_noSF_curved(
							observation_FS, ctf, proj, Ewald_radius,
							data,
							weight,
							num_threads);
					}
					else
					{
						FourierBackprojection::backprojectSlice_noSF_curved_slow(
							observation_FS, ctf, proj, Ewald_radius,
							data,
							weight,
							num_threads);
					}
				}
				else
				{				
					FourierBackprojection::backprojectSlice_noSF(
						observation_FS, ctf, proj,
						data,
						weight,
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
		Image<RFLOAT> vol_xmipp;		
		vol_xmipp().initZeros(s, s, s);
		vol_xmipp().setXmippOrigin();
		
		MultidimArray<RFLOAT> tau2;
		
		backprojector.reconstruct(vol_xmipp(), 10, false, tau2);
		
		data_div_RS.copyDataAndSizeFrom(vol_xmipp);
	}
	else
	{
		data.writeVtk("data_"+tag+".vtk");
		        
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
