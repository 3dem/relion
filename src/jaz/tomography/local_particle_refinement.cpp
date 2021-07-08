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
#define OUTPUT_SCALE 1e7

LocalParticleRefinement::LocalParticleRefinement(
		ParticleIndex particle_id,
		const ParticleSet& particleSet,
		const Tomogram& tomogram,
		const TomoReferenceMap& reference,
		const BufferedImage<float>& freqWeight,
		const BufferedImage<float>& doseWeight,
		const AberrationsCache& aberrationsCache,
		double dose_cutoff,
		int minFrame,
		int maxFrame)
:
	particle_id(particle_id),
	particleSet(particleSet),
	tomogram(tomogram),
	reference(reference),
	freqWeight(freqWeight),
	doseWeight(doseWeight),
	aberrationsCache(aberrationsCache),
	minFrame(minFrame),
	maxFrame(maxFrame)
{
	const int s = reference.getBoxSize();
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = s * pixelSize;
	const int og = particleSet.getOpticsGroup(particle_id);

	if (this->maxFrame < 0) this->maxFrame = fc-1;

	position = particleSet.getPosition(particle_id);

	const std::vector<d3Vector> trajectory = particleSet.getTrajectoryInPixels(
				particle_id, fc, pixelSize);
	
	isVisible = tomogram.determineVisiblity(trajectory, s/2.0);

	observations = BufferedImage<fComplex>(sh,s,fc);

	std::vector<d4Matrix> tomo_to_image;

	TomoExtraction::extractAt3D_Fourier(
			tomogram.stack, s, 1.0, tomogram, trajectory, isVisible,
			observations, tomo_to_image, 1, false);

	const d4Matrix particle_to_tomo = particleSet.getMatrix4x4(
			particle_id, s, s, s);

	Pt.resize(fc);
	CTFs.resize(fc);
	max_radius = std::vector<int>(fc,sh);

	precomputedCTFs = BufferedImage<float>(sh,s,fc);

	for (int f = 0; f < fc; f++)
	{
		if (!isVisible[f]) continue;
		
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

			if (dose_weight < dose_cutoff)
			{
				max_radius[f] = x;
				break;
			}
		}

		if (max_radius[f] > reference.lastShell)
		{
			max_radius[f] = reference.lastShell;
		}

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

					const double gamma_offset = aberrationsCache.hasSymmetrical?
						aberrationsCache.symmetrical[og](xi,yi) : 0.0;

					precomputedCTFs(xi,yi,f) = -doseWeight(xi,yi,f) * CTFs[f].getCTF(
						xA, yA, false, false, false, true, gamma_offset);

				}
			}
		}
	}

	if (aberrationsCache.hasAntisymmetrical)
	{
		aberrationsCache.correctObservations(observations, og);
	}
}

double LocalParticleRefinement::f(const std::vector<double>& x, void* tempStorage) const
{
	for (int i = 0; i < 6; i++)
	{
		if (isnan(x[i])) return std::numeric_limits<double>::max();
	}

	const int s = reference.getBoxSize();
	const int fc = tomogram.frameCount;
	const int hs = particleSet.getHalfSet(particle_id);
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = s * pixelSize;
	const double og = particleSet.getOpticsGroup(particle_id);

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
		if (!isVisible[f]) continue;
		
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

					/*const double gamma_offset = aberrationsCache.hasSymmetrical?
						aberrationsCache.symmetrical[og](xi,yi) : 0.0;

					const float c = -doseWeight(xi,yi,f) * CTFs[f].getCTF(
						xA, yA, false, false, false, true, gamma_offset);*/

					const float c = precomputedCTFs(xi,yi,f);

					const double t = 2.0 * PI * (tx * xp + ty * yp) / (double) s;

					const fComplex shift(cos(t), sin(t));

					const fComplex dF = c * shift * pred - obs;

					L2 += freqWeight(xi,yi,f) * dF.norm();
				}
			}
		}
	}

	const double scale = OUTPUT_SCALE / ((double)s * (double)s * (double)fc);

	return scale * L2;
}

void LocalParticleRefinement::grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	for (int i = 0; i < 6; i++)
	{
		if (isnan(x[i])) return;
	}

	const int s = reference.getBoxSize();
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int hs = particleSet.getHalfSet(particle_id);
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = s * pixelSize;
	const double og = particleSet.getOpticsGroup(particle_id);

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
		if (!isVisible[f]) continue;
		
		const d3Matrix PAt = At * Pt[f];

		const d3Matrix dPAt_dphi   = dAt_dx[0] * Pt[f];
		const d3Matrix dPAt_dtheta = dAt_dx[1] * Pt[f];
		const d3Matrix dPAt_dchi   = dAt_dx[2] * Pt[f];

		const d4Matrix& P = tomogram.projectionMatrices[f];

		const double tx = P(0,0) * tX + P(0,1) * tY + P(0,2) * tZ;
		const double ty = P(1,0) * tX + P(1,1) * tY + P(1,2) * tZ;

		const double dtx_dtX = P(0,0);
		const double dtx_dtY = P(0,1);
		const double dtx_dtZ = P(0,2);

		const double dty_dtX = P(1,0);
		const double dty_dtY = P(1,1);
		const double dty_dtZ = P(1,2);


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

					const t4Vector<fComplex> dPred_dP3D = Interpolation::linearXYZGradientAndValue_FftwHalf_complex(
						reference.image_FS[hs], p3D.x, p3D.y, p3D.z);

					const fComplex pred = dPred_dP3D.w;

					const fComplex dPred_dphi   = (
						dPred_dP3D.x * dP3D_dphi.x   +
						dPred_dP3D.y * dP3D_dphi.y   +
						dPred_dP3D.z * dP3D_dphi.z );

					const fComplex dPred_dTheta = (
						dPred_dP3D.x * dP3D_dtheta.x +
						dPred_dP3D.y * dP3D_dtheta.y +
						dPred_dP3D.z * dP3D_dtheta.z );

					const fComplex dPred_dChi   = (
						dPred_dP3D.x * dP3D_dchi.x   +
						dPred_dP3D.y * dP3D_dchi.y   +
						dPred_dP3D.z * dP3D_dchi.z );


					const fComplex obs = observations(xi,yi,f);

					/*const double gamma_offset = aberrationsCache.hasSymmetrical?
						aberrationsCache.symmetrical[og](xi,yi) : 0.0;

					const float c = -doseWeight(xi,yi,f) * CTFs[f].getCTF(
						xA, yA, false, false, false, true, gamma_offset);*/

					const float c = precomputedCTFs(xi,yi,f);

					const double t = 2.0 * PI * (tx * xp + ty * yp) / (double) s;
					const double dt_dtx =  2.0 * PI * xp / (double) s;
					const double dt_dty =  2.0 * PI * yp / (double) s;

					const fComplex shift(cos(t), sin(t));
					const fComplex dShift_dt(-sin(t), cos(t));

					const fComplex dShift_dtX = (dt_dtx * dtx_dtX + dt_dty * dty_dtX) * dShift_dt;
					const fComplex dShift_dtY = (dt_dtx * dtx_dtY + dt_dty * dty_dtY) * dShift_dt;
					const fComplex dShift_dtZ = (dt_dtx * dtx_dtZ + dt_dty * dty_dtZ) * dShift_dt;


					const fComplex dF = c * shift * pred - obs;

					const fComplex ddF_dPhi   = c * shift * dPred_dphi;
					const fComplex ddF_dTheta = c * shift * dPred_dTheta;
					const fComplex ddF_dChi   = c * shift * dPred_dChi;

					const fComplex ddF_dShift = c * pred;

					const fComplex ddF_dtX = ddF_dShift * dShift_dtX;
					const fComplex ddF_dtY = ddF_dShift * dShift_dtY;
					const fComplex ddF_dtZ = ddF_dShift * dShift_dtZ;


					// L2 += freqWeight(xi,yi,f) * dF.norm();

					const fComplex dL2_ddF = 2.f * freqWeight(xi,yi,f) * dF;


					gradDest[0] += ANGLE_SCALE *
						(dL2_ddF.real * ddF_dPhi.real   + dL2_ddF.imag * ddF_dPhi.imag);

					gradDest[1] += ANGLE_SCALE *
						(dL2_ddF.real * ddF_dTheta.real + dL2_ddF.imag * ddF_dTheta.imag);

					gradDest[2] += ANGLE_SCALE *
						(dL2_ddF.real * ddF_dChi.real   + dL2_ddF.imag * ddF_dChi.imag);


					gradDest[3] += SHIFT_SCALE *
						(dL2_ddF.real * ddF_dtX.real + dL2_ddF.imag * ddF_dtX.imag);

					gradDest[4] += SHIFT_SCALE *
						(dL2_ddF.real * ddF_dtY.real + dL2_ddF.imag * ddF_dtY.imag);

					gradDest[5] += SHIFT_SCALE *
						(dL2_ddF.real * ddF_dtZ.real + dL2_ddF.imag * ddF_dtZ.imag);
				}
			}
		}
	}

	const double scale = OUTPUT_SCALE / ((double)s * (double)s * (double)fc);

	for (int i = 0; i < 6; i++)
	{
		gradDest[i] *= scale;
	}
}

double LocalParticleRefinement::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
{
	for (int i = 0; i < 6; i++)
	{
		gradDest[i] = 0.0;
	}

	for (int i = 0; i < 6; i++)
	{
		if (isnan(x[i])) return std::numeric_limits<double>::max();
	}

	const int s = reference.getBoxSize();
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int hs = particleSet.getHalfSet(particle_id);
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = s * pixelSize;
	const double og = particleSet.getOpticsGroup(particle_id);

	const double phi   = ANGLE_SCALE * x[0];
	const double theta = ANGLE_SCALE * x[1];
	const double chi   = ANGLE_SCALE * x[2];

	const double tX = SHIFT_SCALE * x[3];
	const double tY = SHIFT_SCALE * x[4];
	const double tZ = SHIFT_SCALE * x[5];

	const t4Vector<d3Matrix> dAt_dx = TaitBryan::anglesToMatrixAndDerivatives(phi, theta, chi);
	const d3Matrix& At = dAt_dx[3];

	double L2 = 0.0;

	for (int f = minFrame; f <= maxFrame; f++)
	{
		if (!isVisible[f]) continue;
		
		const d3Matrix PAt = At * Pt[f];

		const d3Matrix dPAt_dphi   = dAt_dx[0] * Pt[f];
		const d3Matrix dPAt_dtheta = dAt_dx[1] * Pt[f];
		const d3Matrix dPAt_dchi   = dAt_dx[2] * Pt[f];

		const d4Matrix& P = tomogram.projectionMatrices[f];

		const double tx = P(0,0) * tX + P(0,1) * tY + P(0,2) * tZ;
		const double ty = P(1,0) * tX + P(1,1) * tY + P(1,2) * tZ;

		const double dtx_dtX = P(0,0);
		const double dtx_dtY = P(0,1);
		const double dtx_dtZ = P(0,2);

		const double dty_dtX = P(1,0);
		const double dty_dtY = P(1,1);
		const double dty_dtZ = P(1,2);


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

					const t4Vector<fComplex> dPred_dP3D = Interpolation::linearXYZGradientAndValue_FftwHalf_complex(
						reference.image_FS[hs], p3D.x, p3D.y, p3D.z);

					const fComplex pred = dPred_dP3D.w;

					const fComplex dPred_dPhi   = (
						dPred_dP3D.x * dP3D_dphi.x   +
						dPred_dP3D.y * dP3D_dphi.y   +
						dPred_dP3D.z * dP3D_dphi.z );

					const fComplex dPred_dTheta = (
						dPred_dP3D.x * dP3D_dtheta.x +
						dPred_dP3D.y * dP3D_dtheta.y +
						dPred_dP3D.z * dP3D_dtheta.z );

					const fComplex dPred_dChi   = (
						dPred_dP3D.x * dP3D_dchi.x   +
						dPred_dP3D.y * dP3D_dchi.y   +
						dPred_dP3D.z * dP3D_dchi.z );


					const fComplex obs = observations(xi,yi,f);

					/*const double gamma_offset = aberrationsCache.hasSymmetrical?
						aberrationsCache.symmetrical[og](xi,yi) : 0.0;

					const float c = -doseWeight(xi,yi,f) * CTFs[f].getCTF(
						xA, yA, false, false, false, true, gamma_offset);*/

					const float c = precomputedCTFs(xi,yi,f);

					const double t = 2.0 * PI * (tx * xp + ty * yp) / (double) s;
					const double dt_dtx =  2.0 * PI * xp / (double) s;
					const double dt_dty =  2.0 * PI * yp / (double) s;

					const fComplex shift(cos(t), sin(t));
					const fComplex dShift_dt(-sin(t), cos(t));

					const fComplex dShift_dtX = (dt_dtx * dtx_dtX + dt_dty * dty_dtX) * dShift_dt;
					const fComplex dShift_dtY = (dt_dtx * dtx_dtY + dt_dty * dty_dtY) * dShift_dt;
					const fComplex dShift_dtZ = (dt_dtx * dtx_dtZ + dt_dty * dty_dtZ) * dShift_dt;


					const fComplex dF = c * shift * pred - obs;

					const fComplex ddF_dPhi   = c * shift * dPred_dPhi;
					const fComplex ddF_dTheta = c * shift * dPred_dTheta;
					const fComplex ddF_dChi   = c * shift * dPred_dChi;

					const fComplex ddF_dShift = c * pred;

					const fComplex ddF_dtX = ddF_dShift * dShift_dtX;
					const fComplex ddF_dtY = ddF_dShift * dShift_dtY;
					const fComplex ddF_dtZ = ddF_dShift * dShift_dtZ;


					L2 += freqWeight(xi,yi,f) * dF.norm();

					const fComplex dL2_ddF = 2.f * freqWeight(xi,yi,f) * dF;


					gradDest[0] += ANGLE_SCALE *
						(dL2_ddF.real * ddF_dPhi.real   + dL2_ddF.imag * ddF_dPhi.imag);

					gradDest[1] += ANGLE_SCALE *
						(dL2_ddF.real * ddF_dTheta.real + dL2_ddF.imag * ddF_dTheta.imag);

					gradDest[2] += ANGLE_SCALE *
						(dL2_ddF.real * ddF_dChi.real   + dL2_ddF.imag * ddF_dChi.imag);


					gradDest[3] += SHIFT_SCALE *
						(dL2_ddF.real * ddF_dtX.real + dL2_ddF.imag * ddF_dtX.imag);

					gradDest[4] += SHIFT_SCALE *
						(dL2_ddF.real * ddF_dtY.real + dL2_ddF.imag * ddF_dtY.imag);

					gradDest[5] += SHIFT_SCALE *
						(dL2_ddF.real * ddF_dtZ.real + dL2_ddF.imag * ddF_dtZ.imag);
				}
			}
		}
	}

	const double scale = OUTPUT_SCALE / ((double)s * (double)s * (double)fc);

	for (int i = 0; i < 6; i++)
	{
		gradDest[i] *= scale;
	}

	return scale * L2;
}

void LocalParticleRefinement::applyChange(const std::vector<double>& x, ParticleSet& target, ParticleIndex particle_id, double pixel_size)
{
	d3Matrix A0 = target.getParticleMatrix(particle_id);
	d3Matrix B = TaitBryan::anglesToMatrix3(ANGLE_SCALE * x[0], ANGLE_SCALE * x[1], ANGLE_SCALE * x[2]).transpose();
	d3Matrix A = A0 * B;

	d3Vector ang = Euler::matrixToAngles(A);

	target.partTable.setValue(EMDL_ORIENT_ROT, RAD2DEG(ang[0]), particle_id.value);
	target.partTable.setValue(EMDL_ORIENT_TILT, RAD2DEG(ang[1]), particle_id.value);
	target.partTable.setValue(EMDL_ORIENT_PSI, RAD2DEG(ang[2]), particle_id.value);

	d3Vector pos0;

	pos0[0] = target.partTable.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, particle_id.value);
	pos0[1] = target.partTable.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, particle_id.value);
	pos0[2] = target.partTable.getDouble(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, particle_id.value);

	target.partTable.setValue(
			EMDL_ORIENT_ORIGIN_X_ANGSTROM,
			pos0[0] + pixel_size * SHIFT_SCALE * x[3],
			particle_id.value);

	target.partTable.setValue(
			EMDL_ORIENT_ORIGIN_Y_ANGSTROM,
			pos0[1] + pixel_size * SHIFT_SCALE * x[4],
			particle_id.value);

	target.partTable.setValue(
			EMDL_ORIENT_ORIGIN_Z_ANGSTROM,
			pos0[2] + pixel_size * SHIFT_SCALE * x[5],
			particle_id.value);
}
