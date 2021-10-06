#include "tomo_mag_fit.h"
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/projection/fwd_projection.h>
#include <omp.h>

using namespace gravis;


TomoMagFit::TomoMagFit(
		const std::vector<ParticleIndex>& particle_indices,
		const Tomogram& tomogram,
		const ParticleSet& particleSet,
		const TomoReferenceMap& referenceMap,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		int boxSize,
		int first_frame,
		int last_frame,
		int num_threads)
:
	particle_indices(particle_indices),
	tomogram(tomogram),
	particleSet(particleSet),
	referenceMap(referenceMap),
	freqWeights(freqWeights),
	doseWeights(doseWeights),
	boxSize(boxSize),
	first_frame(first_frame),
	last_frame(last_frame),
	num_threads(num_threads)
{
}

TomoIsoMagFit::TomoIsoMagFit(
		const std::vector<ParticleIndex>& particle_indices,
		const Tomogram& tomogram,
		const ParticleSet& particleSet,
		const TomoReferenceMap& referenceMap,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		int boxSize,
		int first_frame,
		int last_frame,
		int num_threads)
:
	TomoMagFit(
		particle_indices, tomogram, particleSet, referenceMap,
		freqWeights, doseWeights,
		boxSize, first_frame, last_frame, num_threads)
{
}

d2Vector TomoIsoMagFit::computeErrorAndSlope(
		double mag,
		bool consider_image_scale,
		bool consider_defocus_stretch)
{
	const int pc = particle_indices.size();
	const int fc = tomogram.frameCount;

	const double pixelSize0 = tomogram.optics.pixelSize;
	const double ba0 = pixelSize0 * boxSize;

	const int s = boxSize;
	const int sh = s/2 + 1;

	d3Vector centre_of_mass(0.0, 0.0, 0.0);

	for (int p = 0; p < pc; p++)
	{
		const ParticleIndex particle_id = particle_indices[p];
		const d3Vector pos = particleSet.getPosition(particle_id);
		centre_of_mass += pos;
	}

	centre_of_mass /= pc;

	std::vector<double> mean_depth(fc, 0.0);

	for (int f = first_frame; f <= last_frame; f++)
	{
		mean_depth[f] = tomogram.getDepthOffset(f, centre_of_mass);
	}


	const int data_pad = 256;
	std::vector<d2Vector> out_per_thread(num_threads * data_pad, d2Vector(0.0, 0.0));

	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		const ParticleIndex particle_id = particle_indices[p];

		const d3Vector pos = particleSet.getPosition(particle_id);
		const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
					particle_id, fc, pixelSize0);
		
		const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);

		d4Matrix projCut;

		BufferedImage<fComplex> observation(sh,s);

		double L2 = 0.0;
		double dL2_dmag = 0.0;


		for (int f = first_frame; f <= last_frame; f++)
		{
			if (!isVisible[f]) continue;
			
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram,
					traj[f], observation, projCut, 1, true);


			const d4Matrix particleToTomo = particleSet.getMatrix4x4(
					particle_id, s, s, s);

			d4Matrix projPart;

			if (consider_image_scale)
			{
				projPart = mag * projCut * particleToTomo;
			}
			else
			{
				projPart = projCut * particleToTomo;
			}

			const int hs = particleSet.getHalfSet(particle_id);

			BufferedImage<fComplex> prediction(sh,s);
			BufferedImage<t2Vector<fComplex>> predGradient(sh,s);

			ForwardProjection::forwardProject(
					referenceMap.image_FS[hs], {projPart}, prediction, 1);

			ForwardProjection::forwardProject2DGradient(
					referenceMap.image_FS[hs], {projPart}, predGradient, 1);

			const float scale = -1.f;

			// construct ctf with correct stretching factor
			// (ensure that the defocus of the centre of mass does not change)

			const double dz0 = tomogram.getDepthOffset(f, pos);

			double dz, ddeltaF_dinvmag;

			if (consider_defocus_stretch)
			{
				dz  = mean_depth[f] + (dz0 - mean_depth[f]) / mag;
				ddeltaF_dinvmag = tomogram.handedness * tomogram.optics.pixelSize * (dz0 - mean_depth[f]);
			}
			else
			{
				dz = dz0;
				ddeltaF_dinvmag = 0.0;
			}

			const double deltaF = tomogram.handedness * tomogram.optics.pixelSize * dz;
			const double ddeltaF_dmag = -ddeltaF_dinvmag / (mag * mag);

			CTF ctf_part = tomogram.centralCTFs[f];

			ctf_part.DeltafU += deltaF;
			ctf_part.DeltafV += deltaF;

			ctf_part.initialise();

			std::vector<double> K_ctf = ctf_part.getK();


			// d_gamma / d_DeltaF = -lambda * PI * r² = -K_ctf[1] * r²

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const double xp = x;
				const double yp = y < s/2? y : y - s;

				const double xa = xp / ba0;
				const double ya = yp / ba0;

				const double gamma = ctf_part.getLowOrderGamma(xa,ya);

				const float c = -scale * sin(gamma) * doseWeights(x,y,f);

				const float wg = freqWeights(x,y,f);

				const fComplex pred = ctf_part.scale * prediction(x,y);
				const t2Vector<fComplex> grad = predGradient(x,y);

				const fComplex dF = c * pred - observation(x,y);

				const float dgamma_ddeltaF = -K_ctf[1] * (xa*xa + ya*ya);
				const float dgamma_dmag = dgamma_ddeltaF * ddeltaF_dmag;

				const float dc_dmag = -scale * cos(gamma) * dgamma_dmag;

				const fComplex dpred_dmag = fComplex(
						grad.x.real * xp + grad.y.real * yp,
						grad.x.imag * xp + grad.y.imag * yp);

				fComplex dFp_dmag(0.f, 0.f);

				if (consider_image_scale)
				{
					dFp_dmag += c * dpred_dmag;
				}

				if (consider_defocus_stretch)
				{
					dFp_dmag += dc_dmag * pred;
				}

				L2 += wg * dF.norm();

				dL2_dmag += 2.0 * wg * (dF.real * dFp_dmag.real + dF.imag * dFp_dmag.imag);
			}

			out_per_thread[data_pad * th] += d2Vector(L2, dL2_dmag);
		}
	}

	d2Vector out(0.0, 0.0);

	for (int th = 0; th < num_threads; th++)
	{
		out += out_per_thread[data_pad * th];
	}

	return out;
}

TomoAnisoMagFit::TomoAnisoMagFit(
		const std::vector<ParticleIndex>& particle_indices,
		const Tomogram& tomogram,
		const ParticleSet& particleSet,
		const TomoReferenceMap& referenceMap,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		int boxSize,
		int first_frame,
		int last_frame,
		int num_threads)
:
	TomoMagFit(
		particle_indices, tomogram, particleSet, referenceMap,
		freqWeights, doseWeights,
		boxSize, first_frame, last_frame, num_threads)
{

}

BufferedImage<Equation2x2> TomoAnisoMagFit::computeEquations()
{
	const int pc = particle_indices.size();
	const int fc = tomogram.frameCount;

	const double pixelSize = tomogram.optics.pixelSize;

	const int s = boxSize;
	const int sh = s/2 + 1;

	const int data_pad = 2;

	std::vector<BufferedImage<Equation2x2>> equations_per_thread(
					data_pad * num_threads, BufferedImage<Equation2x2>(sh,s));


	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		const ParticleIndex particle_id = particle_indices[p];

		const d3Vector pos = particleSet.getPosition(particle_id);
		const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
					particle_id, fc, pixelSize);
		
		const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);

		d4Matrix projCut;

		BufferedImage<fComplex> observation(sh,s);

		for (int f = first_frame; f <= last_frame; f++)
		{
			if (!isVisible[f]) continue;
			
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram,
					traj[f], observation, projCut, 1, true);

			observation *= -1.f;

			const d4Matrix particleToTomo = particleSet.getMatrix4x4(
					particle_id, s, s, s);

			d4Matrix projPart = projCut * particleToTomo;

			const int hs = particleSet.getHalfSet(particle_id);

			BufferedImage<fComplex> prediction(sh,s);
			BufferedImage<t2Vector<fComplex>> predGradient(sh,s);

			ForwardProjection::forwardProject(
					referenceMap.image_FS[hs], {projPart}, prediction, 1);

			ForwardProjection::forwardProject2DGradient(
					referenceMap.image_FS[hs], {projPart}, predGradient, 1);


			const double dz0 = tomogram.getDepthOffset(f, pos);
			const double deltaF = tomogram.handedness * tomogram.optics.pixelSize * dz0;

			CTF ctf_part = tomogram.centralCTFs[f];

			ctf_part.DeltafU += deltaF;
			ctf_part.DeltafV += deltaF;

			ctf_part.initialise();


			MagnificationHelper::updateScale(
					prediction, predGradient, observation,
					freqWeights.getConstSliceRef(f),
					doseWeights.getConstSliceRef(f),
					ctf_part, pixelSize, equations_per_thread[data_pad * th]);
		}
	}

	BufferedImage<Equation2x2> equations(sh,s);

	for (int th = 0; th < num_threads; th++)
	{
		equations += equations_per_thread[data_pad * th];
	}

	return equations;
}

std::vector<BufferedImage<Equation2x2>> TomoAnisoMagFit::computeEquations_even_odd()
{
	const int pc = particle_indices.size();
	const int fc = tomogram.frameCount;

	const double pixelSize = tomogram.optics.pixelSize;

	const int s = boxSize;
	const int sh = s/2 + 1;

	const int data_pad = 2;

	std::vector<std::vector<BufferedImage<Equation2x2>>> equations_per_thread(
				3, std::vector<BufferedImage<Equation2x2>>(
					data_pad * num_threads, BufferedImage<Equation2x2>(sh,s)));

	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		const ParticleIndex particle_id = particle_indices[p];

		const d3Vector pos = particleSet.getPosition(particle_id);
		const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
					particle_id, fc, pixelSize);
		
		const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);

		d4Matrix projCut;

		BufferedImage<fComplex> observation(sh,s);

		for (int f = first_frame; f <= last_frame; f++)
		{
			if (!isVisible[f]) continue;
			
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram,
					traj[f], observation, projCut, 1, true);

			observation *= -1.f;

			const d4Matrix particleToTomo = particleSet.getMatrix4x4(
					particle_id, s, s, s);

			d4Matrix projPart = projCut * particleToTomo;

			const int hs = particleSet.getHalfSet(particle_id);

			BufferedImage<fComplex> prediction(sh,s);
			BufferedImage<t2Vector<fComplex>> predGradient(sh,s);

			ForwardProjection::forwardProject(
					referenceMap.image_FS[hs], {projPart}, prediction, 1);

			ForwardProjection::forwardProject2DGradient(
					referenceMap.image_FS[hs], {projPart}, predGradient, 1);


			const double dz0 = tomogram.getDepthOffset(f, pos);
			const double deltaF = tomogram.handedness * tomogram.optics.pixelSize * dz0;

			CTF ctf_part = tomogram.centralCTFs[f];

			ctf_part.DeltafU += deltaF;
			ctf_part.DeltafV += deltaF;

			ctf_part.initialise();


			MagnificationHelper::updateScale(
					prediction, predGradient, observation,
					freqWeights.getConstSliceRef(f),
					doseWeights.getConstSliceRef(f),
					ctf_part, pixelSize, equations_per_thread[0][data_pad * th]);

			if (f%2 == 0)
			{
				MagnificationHelper::updateScale(
						prediction, predGradient, observation,
						freqWeights.getConstSliceRef(f),
						doseWeights.getConstSliceRef(f),
						ctf_part, pixelSize, equations_per_thread[1][data_pad * th]);
			}
			else
			{
				MagnificationHelper::updateScale(
						prediction, predGradient, observation,
						freqWeights.getConstSliceRef(f),
						doseWeights.getConstSliceRef(f),
						ctf_part, pixelSize, equations_per_thread[2][data_pad * th]);
			}
		}
	}

	std::vector<BufferedImage<Equation2x2>> equations(3,
				BufferedImage<Equation2x2>(sh,s));

	for (int i = 0; i < 3; i++)
	{
		for (int th = 0; th < num_threads; th++)
		{
			equations[i] += equations_per_thread[i][data_pad * th];
		}
	}

	return equations;
}

BufferedImage<double> TomoAnisoMagFit::computePerPixelSlope(const RawImage<Equation2x2> &equations)
{
	return MagnificationHelper::solvePerPixel(equations);
}

double TomoAnisoMagFit::evaluateMag(const d2Matrix& M)
{
	const int pc = particle_indices.size();
	const int fc = tomogram.frameCount;

	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = pixelSize * boxSize;

	const int s = boxSize;
	const int sh = s/2 + 1;

	const int data_pad = 256;
	const float scale = -1;

	std::vector<double> L2_per_thread(num_threads * data_pad, 0.0);

	const d4Matrix M4t(
			M(0,0), M(1,0), 0.0, 0.0,
			M(0,1), M(1,1), 0.0, 0.0,
			   0.0,    0.0, 1.0, 0.0,
			   0.0,    0.0, 0.0, 1.0);

	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		const ParticleIndex particle_id = particle_indices[p];

		const d3Vector pos = particleSet.getPosition(particle_id);
		const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
					particle_id, fc, pixelSize);
		
		const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);

		d4Matrix projCut;

		BufferedImage<fComplex> observation(sh,s);

		for (int f = first_frame; f <= last_frame; f++)
		{
			if (!isVisible[f]) continue;
			
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram,
					traj[f], observation, projCut, 1, true);

			const d4Matrix particleToTomo = particleSet.getMatrix4x4(
					particle_id, s, s, s);

			d4Matrix projPart = M4t * projCut * particleToTomo;

			const int hs = particleSet.getHalfSet(particle_id);

			BufferedImage<fComplex> prediction(sh,s);

			ForwardProjection::forwardProject(
					referenceMap.image_FS[hs], {projPart}, prediction, 1);

			const double dz0 = tomogram.getDepthOffset(f, pos);
			const double deltaF = tomogram.handedness * tomogram.optics.pixelSize * dz0;

			CTF ctf_part = tomogram.centralCTFs[f];

			ctf_part.DeltafU += deltaF;
			ctf_part.DeltafV += deltaF;

			ctf_part.initialise();

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const double xp = x;
				const double yp = y < s/2? y : y - s;

				const double xa = xp / ba;
				const double ya = yp / ba;

				const double gamma = ctf_part.getLowOrderGamma(xa,ya);

				const float c = -scale * sin(gamma) * doseWeights(x,y,f);

				const float wg = freqWeights(x,y,f);

				const fComplex pred = prediction(x,y);
				const fComplex dF = c * pred - observation(x,y);

				L2_per_thread[th * data_pad] += wg * dF.norm();
			}
		}
	}

	double L2 = 0.0;

	for (int th = 0; th < num_threads; th++)
	{
		L2 += L2_per_thread[data_pad * th];
	}

	return L2;
}

