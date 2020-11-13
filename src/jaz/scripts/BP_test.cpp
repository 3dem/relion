#include <src/args.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/projection/point_insertion.h>
#include <src/jaz/tomography/projection/fwd_projection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/reference_map.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/math/Euler_angles_dynamo.h>
#include <src/CPlot2D.h>
#include <src/backprojector.h>

using namespace gravis;

#define INPUT_PRECISION float
#define OUTPUT_PRECISION float

int main(int argc, char *argv[])
{
	const int s = 180;
	const int sh = s/2 + 1;
	const int num_observations = 100000;
	const int num_shells = 5;
	const int num_spheres = 2 * num_shells - 1;
	const int num_threads = 6;
	const double outer_radius = s / 2;
	
	const std::string tag = "100k_backward_noxg";
	
	const double SNR = 1.0;
	const bool forward = false;
	const bool explicit_gridding = false;
	const bool legacy_backprojector = false;
	const bool wrap_voxels = false;
	
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
		
		d4Matrix proj = EulerDynamo::anglesToMatrix4(phi, tilt, psi);
				
		if (legacy_backprojector)
		{
			Matrix2D<RFLOAT> A(3,3);
			
			for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++)
			{
				A(r,c) = proj(r,c);
			}
			
			backprojector.set2DFourierTransform(obervation_FS_legacy(), A, &ctf_legacy());
		}
		else
		{
			if (forward)
			{
				if (wrap_voxels)
				{
					FourierBackprojection::backprojectSlice_forward_wrap(
						observation_FS, ctf, proj,
						data, weight);
				}
				else
				{
					ClippedPointInsertion<INPUT_PRECISION,OUTPUT_PRECISION> clippedInsertion;
					
					FourierBackprojection::backprojectSlice_forward(
						clippedInsertion, observation_FS, ctf, proj,
						data, weight);
				}
			}
			else
			{
				FourierBackprojection::backprojectSlice_backward(
					observation_FS, ctf, proj,
					data, weight, num_threads);
				
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
			std::cout << "Reconstruction::griddingCorrect3D_sinc2" << std::endl;
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

	Reconstruction::taper(data_div_RS, 10, true, num_threads);
	
	data_div_RS.write("reconstruction_"+tag+".mrc");


	const int num_predictions = 10;

	BufferedImage<float> predictions(s,s,num_predictions);

	for (int i = 0; i < num_predictions; i++)
	{
		BufferedImage<tComplex<OUTPUT_PRECISION>> data_div_FS;
		FFT::FourierTransform(data_div_RS, data_div_FS);

		Centering::shiftInSitu(data_div_FS);

		const double phi = 2.0 * PI * rand() / (double)RAND_MAX;
		const double psi = 2.0 * PI * rand() / (double)RAND_MAX;
		const double tilt = (PI/2.0) * sin(PI * rand() / (double)RAND_MAX - PI/2.0);

		d4Matrix proj = EulerDynamo::anglesToMatrix4(phi, tilt, psi);


		BufferedImage<fComplex> prediction(sh,s), psf(sh,s);

		ForwardProjection::forwardProject_withPSF(data_div_FS, {proj}, prediction, psf, 1);

		BufferedImage<float> predictionReal(s,s);

		Reconstruction::correctStack(prediction, psf, predictionReal, true, 1);

		predictions.getSliceRef(i).copyFrom(predictionReal);
	}

	predictions.write("predictions_"+tag+".mrc");
	
	std::vector<double> radial_mean(sh, 0.0);
	std::vector<double> radial_count(sh, 0.0);
	
	for (long int z = 0; z < s; z++)
	for (long int y = 0; y < s; y++)
	for (long int x = 0; x < s; x++)
	{
		const double xx = x - s/2;
		const double yy = y - s/2;
		const double zz = z - s/2;
		
		const double r = sqrt(xx*xx + yy*yy + zz*zz);
		const int ri = (int) (r + 0.5);
		
		if (ri < sh)
		{
			radial_mean[ri] += data_div_RS(x,y,z);
			radial_count[ri] += 1;
		}
	}
	
	for (int r = 0; r < sh; r++)
	{
		if (radial_count[r] > 0.0)
		{
			radial_mean[r] /= radial_count[r];
		}
	}
	
	std::vector<double> radial_variance(sh, 0.0);
	
	for (long int z = 0; z < s; z++)
	for (long int y = 0; y < s; y++)
	for (long int x = 0; x < s; x++)
	{
		const double xx = x - s/2;
		const double yy = y - s/2;
		const double zz = z - s/2;
		
		const double r = sqrt(xx*xx + yy*yy + zz*zz);
		const int ri = (int) (r + 0.5);
		
		if (ri < sh)
		{
			const double d = data_div_RS(x,y,z) - radial_mean[r];
			
			radial_variance[ri] += d * d;
		}
	}
	
	for (int r = 0; r < sh; r++)
	{
		if (radial_count[r] > 1.0)
		{
			radial_variance[r] /= radial_count[r] - 1.0;
		}
	}
	
	double plot_mean = 0.0;
	
	for (int r = 0; r < sh; r++)
	{
		plot_mean += radial_mean[r];
	}
	
	plot_mean /= sh;
	
	double plot_power = 0.0;
	
	for (int r = 0; r < sh; r++)
	{
		const double d = radial_mean[r] - plot_mean;
		plot_power += d * d;
	}
	
	const double plot_scale = 1.0 / sqrt(plot_power/sh);
	
		
	std::vector<CDataSet> curves(3);

	for (int i = 0; i < 3; i++)
	{
		CDataSet& curve = curves[i];
		
		curve.SetDrawMarker(false);
		curve.SetDrawLine(true);
		curve.SetLineWidth(0.5);
		
	}
	
	curves[0].SetDatasetColor(0,0,0);
	curves[1].SetDatasetColor(0.5,0.5,0.5);
	curves[2].SetDatasetColor(0.5,0.5,0.5);
	
	for (int r = 0; r < sh; r++)
	{
		const double mu = plot_scale * (radial_mean[r] - plot_mean);
		const double sd = plot_scale * sqrt(radial_variance[r]);
		
		curves[0].AddDataPoint(CDataPoint(r, mu));
		curves[1].AddDataPoint(CDataPoint(r, mu - sd));
		curves[2].AddDataPoint(CDataPoint(r, mu + sd));
	}
	
	CPlot2D plot(tag);
	
	for (int i = 2; i >= 0; i--)
	{
		plot.AddDataSet(curves[i]);
	}
		
	std::string title = tag;
	
	for (int i = 0; i < title.length(); i++)
	{
		if (title[i] == '_')
		{
			title[i] = ' ';
		}
	}
	
	plot.SetTitle(title);
	plot.SetXAxisTitle("r");
	plot.SetYAxisTitle("mean");
	
	plot.SetViewArea(0,-2,sh,3);

	plot.OutputPostScriptPlot("reconstruction_"+tag+".eps");
	
	return 0;
}
