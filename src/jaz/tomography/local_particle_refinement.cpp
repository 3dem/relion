#include "local_particle_refinement.h"
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/projection/fwd_projection.h>
#include <src/jaz/math/Tait_Bryan_angles.h>
#include <src/jaz/math/Euler_angles_relion.h>
#include <src/jaz/optics/damage.h>

using namespace gravis;

#define ANGLE_SCALE 0.02
#define SHIFT_SCALE 1.0

LocalParticleRefinement::LocalParticleRefinement(
		ParticleIndex particle_id,
		const ParticleSet& particleSet,
		const Tomogram& tomogram,
		const TomoReferenceMap& reference,
		const BufferedImage<float>& frqWeight,
		const AberrationsCache& aberrationsCache,
		bool debug)
:
	particle_id(particle_id),
	particleSet(particleSet),
	tomogram(tomogram),
	reference(reference),
	frqWeight(frqWeight),
	aberrationsCache(aberrationsCache),
	debug(debug)
{
	const int s = reference.getBoxSize();
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = s * pixelSize;

	position = particleSet.getPosition(particle_id);

	const std::vector<d3Vector> trajectory = particleSet.getTrajectoryInPixels(
				particle_id, fc, pixelSize);

	observations = BufferedImage<fComplex>(sh,s,fc);

	std::vector<d4Matrix> tomo_to_image;

	TomoExtraction::extractAt3D_Fourier(
			tomogram.stack, s, 1.0, tomogram.projectionMatrices,
				trajectory, observations, tomo_to_image, 1, true);

	const d4Matrix particle_to_tomo = particleSet.getMatrix4x4(
			particle_id, s, s, s);

	Pt.resize(fc);
	CTFs.resize(fc);
	max_radius = std::vector<int>(fc,sh);

	const double weight_threshold = 0.05;

	for (int f = 0; f < fc; f++)
	{
		const d4Matrix A = tomo_to_image[f] * particle_to_tomo;

		Pt[f] = d3Matrix(
					A(0,0), A(1,0), A(2,0),
					A(0,1), A(1,1), A(2,1),
					A(0,2), A(1,2), A(2,2) );

		CTFs[f] = tomogram.getCtf(f, position);

		for (int x = 0; x < sh; x++)
		{
			const double dose = tomogram.cumulativeDose[f];
			const double dose_weight = Damage::getWeight(dose, x / ba);

			if (dose_weight < weight_threshold)
			{
				max_radius[f] = x;
				break;
			}
		}
	}
}

double LocalParticleRefinement::f(const std::vector<double>& x, void* tempStorage) const
{
	const int s = reference.getBoxSize();
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int hs = particleSet.getHalfSet(particle_id);
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = s * pixelSize;

	// @TODO: consider aberrations
	const double gamma_offset = 0.0;

	const double phi   = ANGLE_SCALE * x[0];
	const double theta = ANGLE_SCALE * x[1];
	const double chi   = ANGLE_SCALE * x[2];

	const double tX = SHIFT_SCALE * x[3];
	const double tY = SHIFT_SCALE * x[4];
	const double tZ = SHIFT_SCALE * x[5];

	const d3Matrix At = TaitBryan::anglesToMatrix3(phi, theta, chi);

	double L2 = 0.0;

	for (int f = 0; f < fc; f++)
	{
		const d3Matrix PAt = At * Pt[f];

		const d4Matrix& P = tomogram.projectionMatrices[f];

		const double tx = P(0,0) * tX + P(0,1) * tY + P(0,2) * tZ;
		const double ty = P(1,0) * tX + P(1,1) * tY + P(1,2) * tZ;

		const int rad = max_radius[f];

		for (int yi = 0; yi < s;  yi++)
		{
			const double yp = yi < s/2? yi : yi - s;

			if (yp < rad && yp > -rad)
			{
				const int x_max = (int) sqrt(rad*rad - yp*yp);

				for (int xi = 0; xi < x_max; xi++)
				{
					const double xp = xi;

					const double xA = xp / ba;
					const double yA = yp / ba;

					const d3Vector p2D(xp, yp, 0.0);
					const d3Vector p3D = PAt * p2D;

					const fComplex pred = Interpolation::linearXYZ_FftwHalf_complex(
						reference.image_FS[hs], p3D.x, p3D.y, p3D.z);

					const fComplex obs = observations(xi,yi,f);

					const float c = -CTFs[f].getCTF(
						xA, yA, false, false, false, true, gamma_offset);

					const double t = 2.0 * PI * (tx * xp + ty * yp) / (double) s;

					const fComplex shift(cos(t), sin(t));

					const fComplex dF = c * shift * pred - obs;

					L2 += frqWeight(xi,yi,f) * dF.norm();
				}
			}
		}
	}

	return L2;
}

void LocalParticleRefinement::grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	const int s = reference.getBoxSize();
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int hs = particleSet.getHalfSet(particle_id);
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = s * pixelSize;

	// @TODO: consider aberrations
	const double gamma_offset = 0.0;

	const double phi   = ANGLE_SCALE * x[0];
	const double theta = ANGLE_SCALE * x[1];
	const double chi   = ANGLE_SCALE * x[2];

	const double tX = SHIFT_SCALE * x[3];
	const double tY = SHIFT_SCALE * x[4];
	const double tZ = SHIFT_SCALE * x[5];

	const t4Vector<d3Matrix> dAt_dx = TaitBryan::anglesToMatrixAndDerivatives(phi, theta, chi);
	const d3Matrix& At = dAt_dx[3];

	for (int i = 0; i < 6; i++)
	{
		gradDest[i] = 0.0;
	}


	for (int f = 0; f < fc; f++)
	{
		const d3Matrix PAt = At * Pt[f];

		const d3Matrix dPAt_dphi   = dAt_dx[0] * Pt[f];
		const d3Matrix dPAt_dtheta = dAt_dx[1] * Pt[f];
		const d3Matrix dPAt_dchi   = dAt_dx[2] * Pt[f];

		const d4Matrix& P = tomogram.projectionMatrices[f];

		const double tx = P(0,0) * tX + P(0,1) * tY + P(0,2) * tZ;
		const double ty = P(1,0) * tX + P(1,1) * tY + P(1,2) * tZ;

		const int rad = max_radius[f];

		for (int yi = 0; yi < s;  yi++)
		{
			const double yp = yi < s/2? yi : yi - s;

			if (yp < rad && yp > -rad)
			{
				const int x_max = (int) sqrt(rad*rad - yp*yp);

				for (int xi = 0; xi < x_max; xi++)
				{
					const double xp = xi;

					const double xA = xp / ba;
					const double yA = yp / ba;

					const d3Vector p2D(xp, yp, 0.0);
					const d3Vector p3D = PAt * p2D;

					const d3Vector dP3D_dphi   = dPAt_dphi   * p2D;
					const d3Vector dP3D_dtheta = dPAt_dtheta * p2D;
					const d3Vector dP3D_dchi   = dPAt_dchi   * p2D;

					const fComplex pred = Interpolation::linearXYZ_FftwHalf_complex(
						reference.image_FS[hs], p3D.x, p3D.y, p3D.z);

					const t3Vector<fComplex> dPred_dP3D = Interpolation::linearXYZGradient_FftwHalf_complex(
						reference.image_FS[hs], p3D.x, p3D.y, p3D.z);

					const fComplex obs = observations(xi,yi,f);

					const float c = -CTFs[f].getCTF(
						xA, yA, false, false, false, true, gamma_offset);

					const double t = 2.0 * PI * (tx * xp + ty * yp) / (double) s;

					const fComplex shift(cos(t), sin(t));

					const fComplex dShift_dt(-sin(t), cos(t));
					const float dt_dtx =  2.0 * PI * xp / (float) s;
					const float dt_dty =  2.0 * PI * yp / (float) s;

					const float dtx_dtX = P(0,0);
					const float dtx_dtY = P(0,1);
					const float dtx_dtZ = P(0,2);

					const float dty_dtX = P(1,0);
					const float dty_dtY = P(1,1);
					const float dty_dtZ = P(1,2);

					const fComplex ddF_dShift = c * pred;
					const fComplex dShift_dtX = (dt_dtx * dtx_dtX + dt_dty * dty_dtX) * dShift_dt;
					const fComplex dShift_dtY = (dt_dtx * dtx_dtY + dt_dty * dty_dtY) * dShift_dt;
					const fComplex dShift_dtZ = (dt_dtx * dtx_dtZ + dt_dty * dty_dtZ) * dShift_dt;


					const fComplex dF = c * shift * pred - obs;

					const fComplex ddF_dPhi   = c * shift * (
						dPred_dP3D.x * dP3D_dphi.x   + dPred_dP3D.y * dP3D_dphi.y   + dPred_dP3D.z * dP3D_dphi.z);

					const fComplex ddF_dTheta = c * shift * (
						dPred_dP3D.x * dP3D_dtheta.x + dPred_dP3D.y * dP3D_dtheta.y + dPred_dP3D.z * dP3D_dtheta.z);

					const fComplex ddF_dChi   = c * shift * (
						dPred_dP3D.x * dP3D_dchi.x   + dPred_dP3D.y * dP3D_dchi.y   + dPred_dP3D.z * dP3D_dchi.z);


					const fComplex alpha = 2.f * frqWeight(xi,yi,f) * dF;

					gradDest[0] += ANGLE_SCALE * alpha * (dF.real * ddF_dPhi.real   + dF.imag * ddF_dPhi.imag);
					gradDest[1] += ANGLE_SCALE * alpha * (dF.real * ddF_dTheta.real + dF.imag * ddF_dTheta.imag);
					gradDest[2] += ANGLE_SCALE * alpha * (dF.real * ddF_dChi.real   + dF.imag * ddF_dChi.imag);

					gradDest[3] += SHIFT_SCALE * alpha *
							(dF.real * ddF_dShift.real * dShift_dtX.real +
							 dF.imag * ddF_dShift.imag * dShift_dtX.imag);

					gradDest[4] += SHIFT_SCALE * alpha *
							(dF.real * ddF_dShift.real * dShift_dtY.real +
							 dF.imag * ddF_dShift.imag * dShift_dtY.imag);

					gradDest[5] += SHIFT_SCALE * alpha *
							(dF.real * ddF_dShift.real * dShift_dtZ.real +
							 dF.imag * ddF_dShift.imag * dShift_dtZ.imag);
				}
			}
		}
	}
}

void LocalParticleRefinement::applyChange(const std::vector<double>& x, ParticleSet& target)
{
	d3Matrix A0 = particleSet.getParticleMatrix(particle_id);
	d3Matrix B = TaitBryan::anglesToMatrix3(ANGLE_SCALE * x[0], ANGLE_SCALE * x[1], ANGLE_SCALE * x[2]).transpose();
	d3Matrix A = A0 * B;

	d3Vector ang = Euler::matrixToAngles(A);

	target.partTable.setValue(EMDL_ORIENT_ROT, RAD2DEG(ang[0]), particle_id.value);
	target.partTable.setValue(EMDL_ORIENT_TILT, RAD2DEG(ang[1]), particle_id.value);
	target.partTable.setValue(EMDL_ORIENT_PSI, RAD2DEG(ang[2]), particle_id.value);


	d3Vector pos0;

	pos0[0] = particleSet.partTable.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, particle_id.value);
	pos0[1] = particleSet.partTable.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, particle_id.value);
	pos0[2] = particleSet.partTable.getDouble(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, particle_id.value);

	const double pix2ang = tomogram.optics.pixelSize;

	target.partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, pos0[0] + pix2ang * SHIFT_SCALE * x[3], particle_id.value);
	target.partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, pos0[1] + pix2ang * SHIFT_SCALE * x[4], particle_id.value);
	target.partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, pos0[2] + pix2ang * SHIFT_SCALE * x[5], particle_id.value);
}

