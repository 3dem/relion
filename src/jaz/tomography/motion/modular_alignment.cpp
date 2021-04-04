#include "modular_alignment.h"
#include <src/jaz/util/zio.h>
#include <src/jaz/math/Gaussian_process.h>

using namespace gravis;


LinearMotionModel::LinearMotionModel(
		ModularAlignmentMotionParameters motionParameters,
		ModularAlignmentSettings settings,
		std::vector<d3Vector> initialPos,
		double frameDose,
		double pixelSize,
		bool verbose)
{
	const int pc = initialPos.size();

	double sig_vel_px, sig_div_px;

	if (settings.params_scaled_by_dose)
	{
		sig_vel_px = frameDose * motionParameters.sig_vel / pixelSize;
		sig_div_px = frameDose * motionParameters.sig_div / pixelSize;
	}
	else
	{
		sig_vel_px = motionParameters.sig_vel / pixelSize;
		sig_div_px = motionParameters.sig_div / pixelSize;
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

StaticMotionModel::StaticMotionModel()
{

}
