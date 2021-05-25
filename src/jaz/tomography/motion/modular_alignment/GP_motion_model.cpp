#include "GP_motion_model.h"
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/ctf.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/Tait_Bryan_angles.h>
#include <src/jaz/math/Gaussian_process.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <omp.h>


GPMotionModel::GPMotionModel(
	const ParticleSet& dataSet,
	const std::vector<ParticleIndex>& partIndices,
    const Tomogram& tomogram,  
	MotionParameters motionParameters,   
	Settings settings,
	bool verbose)
{	
	pc = partIndices.size();
	initialPos.resize(pc);
	
	for (int p = 0; p < pc; p++)
	{
		initialPos[p] = dataSet.getPosition(partIndices[p]);
	}
	
	const double pixel_size = tomogram.optics.pixelSize;
	const double frame_dose = tomogram.getFrameDose();
	
	double sig_vel_px, sig_div_px;
	
	if (settings.params_scaled_by_dose)
	{
		sig_vel_px = frame_dose * motionParameters.sig_vel / pixel_size;
		sig_div_px = frame_dose * motionParameters.sig_div / pixel_size;
	}
	else
	{
		sig_vel_px = motionParameters.sig_vel / pixel_size;		
		sig_div_px = motionParameters.sig_div / pixel_size;
	}
	
	if (verbose)
	{
		Log::beginSection("Effective motion parameters:");
		
		Log::print("σ_vel = "+ZIO::itoa(sig_vel_px)+" px/frame");
		Log::print("σ_div = "+ZIO::itoa(sig_div_px)+" px/frame");
		
		Log::endSection();
	}
	
	GpKernel* kernel(0);
	
	if (settings.sqExpKernel)
	{
		kernel = new SquareExponentialKernel(sig_vel_px, sig_div_px);
	}
	else
	{
		kernel = new ExponentialKernel(sig_vel_px, sig_div_px);
	}
	
	Matrix2D<double> C = GaussianProcess::computeCovariance(initialPos, kernel);
	
	delete kernel;
	
	if (verbose)
	{
		Log::print("Decomposing covariance matrix");
	}
	
	GaussianProcess::Basis defBasis = GaussianProcess::getBasis(C, settings.maxEDs);
	
	deformationBasis = defBasis.eigenvectors;
	deformationLambda = defBasis.eigenvalues;
	
	bc = deformationLambda.size();
	
	if (verbose)
	{
		if (bc == pc)
		{
			Log::print("Keeping all " + ZIO::itoa(bc) + " eigendeformations");
		}
		else
		{
			Log::print("Keeping " + ZIO::itoa(bc) + " out of " + ZIO::itoa(pc) + " eigendeformations");
		}
	}
}
