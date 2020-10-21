#include "tomo_mag_fit.h"
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/projection/fwd_projection.h>
#include <omp.h>

using namespace gravis;


TomoMagFit::TomoMagFit(
		const std::vector<int>& particle_indices,
		const Tomogram& tomogram,
		const ParticleSet& particleSet,
		const TomoReferenceMap& referenceMap,
		int boxSize,
		int first_frame,
		int last_frame,
		int num_threads)
:
	particle_indices(particle_indices),
	tomogram(tomogram),
	particleSet(particleSet),
	referenceMap(referenceMap),
	boxSize(boxSize),
	first_frame(first_frame),
	last_frame(last_frame),
	num_threads(num_threads)
{

}

TomoIsoMagFit::TomoIsoMagFit(
		const std::vector<int>& particle_indices,
		const Tomogram& tomogram,
		const ParticleSet& particleSet,
		const TomoReferenceMap& referenceMap,
		int boxSize,
		int first_frame,
		int last_frame,
		int num_threads)
:
	TomoMagFit(
		particle_indices, tomogram, particleSet, referenceMap,
		boxSize, first_frame, last_frame, num_threads)
{

}

d2Vector TomoIsoMagFit::computeErrorAndSlope(double mag)
{
	const int pc = particle_indices.size();
	const int fc = tomogram.frameCount;
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = pixelSize * boxSize;
	const int s = boxSize;
	const int sh = s/2 + 1;
	double avg_offset = 0.0;

	std::vector<double> particle_depth(pc);

	/*for (int p = 0; p < pc; p++)
	{
		const d3Vector pos = particleSet.getPosition(p);
		const double dz = tomogram.getDepthOffset(f, pos);

		avg_offset += dz;

		particle_depth[p] = dz;
	}

	avg_offset /= pc;*/

	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		const int particle_id = particle_indices[p];

		const d3Vector pos = particleSet.getPosition(particle_id);
		const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
					particle_id, fc, pixelSize);

		d4Matrix projCut;

		BufferedImage<fComplex> observation(sh,s);

		for (int f = first_frame; f <= last_frame; f++)
		{
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram.projectionMatrices[f],
					traj[f], observation, projCut, 1, false, true);

			/*BufferedImage<fComplex> prediction = Prediction::predictFS(
					part_id, particleSet, projCut, s, referenceMap.image_FS);*/

			const d4Matrix particleToTomo = particleSet.getMatrix4x4(
						particle_id, s, s, s);

			const d4Matrix projPart = projCut * particleToTomo;

			const int hs = particleSet.getHalfSet(particle_id);

			BufferedImage<fComplex> prediction(sh,s);
			BufferedImage<t2Vector<fComplex>> predGradient(sh,s);

			ForwardProjection::forwardProject(
					referenceMap.image_FS[hs], {projPart}, prediction, 1);

			ForwardProjection::forwardProject2DGradient(
					referenceMap.image_FS[hs], {projPart}, predGradient, 1);

			const float scale = -1.f;

			// make ctf with correct stretching factor

			CTF ctf0 = tomogram.centralCTFs[f];

			CTF ctf_part = ctf0;
			ctf_part.initialise();


			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				/*const double xa = x / ba;
				const double ya = (y < s/2? y : y - s) / ba;

				const float c = scale * ctf.;
				const float wg = freqWeights(x,y,f);

				const double xx = x;
				const double yy = y < s/2? y : y - s;
				const double r2 = xx * xx + yy * yy;

				if (r2 < sh2)
				{
					const fComplex zp = c * prediction(x,y);
					const fComplex zo = observation(x,y);

					CCp += wg * (zp - zo).norm();
				}*/
			}
		}
	}
}
