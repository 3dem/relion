/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/single_particle/image_log.h>

#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/math/Zernike_helper.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optics/aberrations_cache.h>


using namespace aberration;
using namespace gravis;


OriginalBasis AberrationFit::fitBasic(
		Image<RFLOAT> phase, Image<RFLOAT> weight, double angpix)
{
	Matrix2D<RFLOAT> A(5,5);
	Matrix1D<RFLOAT> b(5);

	A.initZeros();
	b.initZeros();

	const int sh = phase.data.xdim;
	const int s = phase.data.ydim;

	const double as = angpix * s;

	OriginalBasis basis;
	std::vector<double> vals(5);

	for (int yi = 0; yi < s; yi++)
	for (int xi = 0; xi < sh; xi++)
	{
		const double x = xi/as;
		const double y = yi > sh? (yi - s)/as: yi/as;

		basis.getBasisValues(x, y, &vals[0]);

		const double v = phase(yi,xi);
		const double w = weight(yi,xi);

		for (int r = 0; r < 5; r++)
		{
			b(r) += w * w * vals[r] * v;

			for (int c = 0; c < 5; c++)
			{
				A(r,c) += w * w * vals[r] * vals[c];
			}
		}
	}

	const double tol = 1e-20;
	Matrix1D<RFLOAT> sol(5);
	solve(A, b, sol, tol);

	for (int i = 0; i < 5; i++)
	{
		basis.coefficients[i] = sol(i);
	}

	return basis;
}

Image<RFLOAT> AberrationFit::draw(AberrationBasis *fit, double angpix, int s)
{
	const int sh = s/2 + 1;
	const double as = angpix * s;

	Image<RFLOAT> vis(sh,s);

	std::vector<double> vals(fit->coefficients.size(), 0.0);

	for (int yi = 0; yi < s; yi++)
	for (int xi = 0; xi < sh; xi++)
	{
		const double x = xi/as;
		const double y = yi > sh? (yi - s)/as: yi/as;

		fit->getBasisValues(x, y, &vals[0]);

		double v = 0.0;

		for (int i = 0; i < 5; i++)
		{
			v += fit->coefficients[i] * vals[i];
		}

		vis(yi,xi) = v;
	}

	return vis;
}


void AberrationFit :: considerParticle(
		ParticleIndex part_id,
		const Tomogram& tomogram,
		const TomoReferenceMap& referenceMap,
		const ParticleSet& dataSet,
		const AberrationsCache& aberrationsCache,
		bool flip_value,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		const BufferedImage<int>& xRanges,
		int f0, int f1,
		BufferedImage<EvenData>& even_out,
		BufferedImage<OddData>& odd_out)
{
	const int s = referenceMap.image_real[0].xdim;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int og = dataSet.getOpticsGroup(part_id);
	const double pix2ang = 1.0 / ((double)s * tomogram.optics.pixelSize);

	if (f1 < 0) f1 = fc - 1;


	const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
				part_id, fc, tomogram.optics.pixelSize);
	
	const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);

	d4Matrix projCut;


	BufferedImage<fComplex> observation(sh,s);

	for (int f = f0; f <= f1; f++)
	{
		if (!isVisible[f]) continue;
		
		TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram, traj[f],
				observation, projCut, 1, true);

		CTF ctf = tomogram.getCtf(f, dataSet.getPosition(part_id));
		const RawImage<float> doseSlice = doseWeights.getConstSliceRef(f);

		BufferedImage<fComplex> prediction = Prediction::predictModulated(
				part_id, dataSet, projCut, s,
				ctf,
				tomogram.optics.pixelSize,
				aberrationsCache,
				referenceMap.image_FS,
				Prediction::OwnHalf,
				Prediction::Unmodulated,
				&doseSlice,
				Prediction::CtfScaled,
				&xRanges(0,f));

		const float scale = flip_value? -1.f : 1.f;

		observation(0,0) = fComplex(0.f, 0.f);
		prediction(0,0) = fComplex(0.f, 0.f);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < xRanges(y,f); x++)
		{
			const double x_ang = pix2ang * x;
			const double y_ang = pix2ang * (y < s/2? y : y - s);

			double gamma = ctf.getLowOrderGamma(x_ang, y_ang);

			if (aberrationsCache.hasSymmetrical)
			{
				gamma += aberrationsCache.symmetrical[og](x,y);
			}

			const double cg = cos(gamma);
			const double sg = sin(gamma);
			const double c = -sg;

			fComplex zobs = observation(x,y);
			fComplex zprd = scale * ctf.scale * prediction(x,y);

			if (aberrationsCache.hasAntisymmetrical)
			{
				const fComplex r = aberrationsCache.phaseShift[og](x,y);
				const fComplex z = zobs;

				zobs.real = z.real * r.real + z.imag * r.imag;
				zobs.imag = z.imag * r.real - z.real * r.imag;
			}

			const double zz = zobs.real * zprd.real + zobs.imag * zprd.imag;
			const double zq = zobs.imag * zprd.real - zobs.real * zprd.imag;
			const double nr = zprd.norm();
			const double wg = freqWeights(x,y,f);


			EvenData& ed = even_out(x,y);

			ed.Axx += wg * nr * sg * sg;
			ed.Axy += wg * nr * cg * sg;
			ed.Ayy += wg * nr * cg * cg;

			ed.bx -= wg * zz * sg;
			ed.by -= wg * zz * cg;


			OddData& od = odd_out(x,y);

			od.a += wg * c * c * nr;

			od.b.real += wg * c * zz;
			od.b.imag += wg * c * zq;
		}
	}
}

EvenSolution AberrationFit::solveEven(
		const BufferedImage<EvenData>& data)
{
	const double eps = 1e-30;
	const int s  = data.ydim;
	const int sh = data.xdim;
	const int fc = data.zdim;

	EvenSolution out;

	out.optimum = BufferedImage<dComplex>(sh,s,fc);
	out.phaseShift = BufferedImage<double>(sh,s,fc);
	out.weight = BufferedImage<Tensor2x2<double>>(sh,s,fc);

	for (int f = 0; f < fc; f++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		EvenData d = data(x,y,f);

		d2Vector b(d.bx, d.by);
		d2Matrix A(d.Axx, d.Axy, d.Axy, d.Ayy);

		const double det = A(0,0) * A(1,1) - A(0,1) * A(1,0);
		const double absDet = std::abs(det);

		if (absDet > eps)
		{
			d2Matrix Ai = A;
			Ai.invert();

			const d2Vector opt = Ai * b;

			out.optimum(x,y,f) = dComplex(opt.x, opt.y);
			out.phaseShift(x,y,f) = std::abs(opt.x) > 0.0? atan2(opt.y, opt.x) : 0.0;
			out.weight(x,y,f) = Tensor2x2<double>(d.Axx, d.Axy, d.Ayy);
		}
		else
		{
			out.optimum(x,y,f) = dComplex(0.0, 0.0);
			out.phaseShift(x,y,f) = 0.0;
			out.weight(x,y,f) = Tensor2x2<double>(0.0, 0.0, 0.0);
		}
	}

	return out;
}

std::vector<double> AberrationFit::fitEven(
		const EvenSolution& solution,
		int n_bands,
		const std::vector<double>& initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{
	const d2Matrix mag(1.0, 0.0, 0.0, 1.0);

	const int cc = Zernike::numberOfEvenCoeffs(n_bands);

	if (initialCoeffs.size() != cc)
	{
		REPORT_ERROR_STR(
			"AberrationFit::solveEven: " << initialCoeffs.size() <<
			" initial coefficient provided, but " << cc << " are required.");
	}

	if (writeImages)
	{
		Centering::fftwHalfToHumanFull(solution.phaseShift).write(prefix + "even_phase_per-pixel.mrc");
	}

	BufferedImage<double> nonlinearFit;

	std::vector<double> coeffs = ZernikeHelper::optimiseEvenZernike(
		solution.optimum,
		solution.weight,
		pixelSize,
		mag,
		n_bands,
		initialCoeffs,
		&nonlinearFit);

	if (writeImages)
	{
		Centering::fftwHalfToHumanFull(nonlinearFit).write(prefix + "even_phase_nonlinear-fit.mrc");
	}

	return coeffs;
}

std::vector<double> AberrationFit::solveAndFitEven(
		const BufferedImage<EvenData> &data,
		int n_bands,
		const std::vector<double> &initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{
	EvenSolution solution = solveEven(data);
	return fitEven(solution, n_bands, initialCoeffs, pixelSize, prefix, writeImages);
}

OddSolution AberrationFit::solveOdd(
		const BufferedImage<OddData>& data)
{
	const int s  = data.ydim;
	const int sh = data.xdim;

	OddSolution out;

	out.optimum = BufferedImage<dComplex>(sh,s);
	out.phaseShift = BufferedImage<double> (sh,s);
	out.weight = BufferedImage<double>(sh,s);

	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		OddData d = data(x,y);

		if (d.a > 0.0)
		{
			out.optimum(x,y) = d.b / d.a;
			out.phaseShift(x,y) = d.b.arg();
			out.weight(x,y) = d.a;
		}
		else
		{
			out.optimum(x,y) = dComplex(0.0, 0.0);
			out.phaseShift(x,y) = 0.0;
			out.weight(x,y) = 0.0;
		}
	}

	return out;
}

std::vector<double> AberrationFit::fitOdd(
		const OddSolution& solution,
		int n_bands,
		const std::vector<double>& initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{
	const d2Matrix mag(1.0, 0.0, 0.0, 1.0);

	const int cc = Zernike::numberOfOddCoeffs(n_bands);

	if (initialCoeffs.size() != cc)
	{
		REPORT_ERROR_STR(
			"AberrationFit::solveOdd: " << initialCoeffs.size() <<
			" initial coefficient provided, but " << cc << " are required.");
	}

	if (writeImages)
	{
		Centering::fftwHalfAntisymmetricalToHumanFull(solution.phaseShift).write(
					prefix + "odd_phase_per-pixel.mrc");
	}

	BufferedImage<double> nonlinearFit;

	std::vector<double> coeffs = ZernikeHelper::optimiseOddZernike(
		solution.optimum,
		solution.weight,
		pixelSize,
		mag,
		n_bands,
		initialCoeffs,
		&nonlinearFit);

	if (writeImages)
	{
		Centering::fftwHalfAntisymmetricalToHumanFull(nonlinearFit).write(
					prefix + "odd_phase_nonlinear-fit.mrc");
	}

	return coeffs;
}

std::vector<double> AberrationFit::solveAndFitOdd(
		const BufferedImage<OddData> &data,
		int n_bands,
		const std::vector<double> &initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{
	OddSolution solution = solveOdd(data);
	return fitOdd(solution, n_bands, initialCoeffs, pixelSize, prefix, writeImages);
}


AberrationBasis::AberrationBasis(int dims)
	:   coefficients(dims, 0.0)
{}

void AberrationBasis::offsetCtf(MetaDataTable &mdt, int particle)
{
	// identical to CTF::read() and CTF::initialize():
	double kV, DeltafU, DeltafV, azimuthal_angle, Cs, scale, Q0, phase_shift;

	if (!mdt.getValue(EMDL_CTF_VOLTAGE, kV, particle)) kV = 200;
	if (!mdt.getValue(EMDL_CTF_DEFOCUSU, DeltafU, particle)) DeltafU = 0;
	if (!mdt.getValue(EMDL_CTF_DEFOCUSV, DeltafV, particle)) DeltafV = DeltafU;
	if (!mdt.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, particle)) azimuthal_angle = 0;
	if (!mdt.getValue(EMDL_CTF_CS, Cs, particle)) Cs = 0;
	//if (!mdt.getValue(EMDL_CTF_SCALEFACTOR, scale, particle)) scale = 1;
	if (!mdt.getValue(EMDL_CTF_Q0, Q0, particle)) Q0 = 0;
	//if (!mdt.getValue(EMDL_CTF_PHASESHIFT, phase_shift, particle)) phase_shift = 0;

	//std::cout << DeltafU << ", " << DeltafV << " @ " << azimuthal_angle << "°, " << Cs << ", " << Q0 << "\n";

	double local_Cs = Cs * 1e7;
	double local_kV = kV * 1e3;
	double rad_azimuth = DEG2RAD(azimuthal_angle);

	double defocus_average   = -(DeltafU + DeltafV) * 0.5;
	double defocus_deviation = -(DeltafU - DeltafV) * 0.5;
	double lambda=12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));

	double K1 = (PI / 2) * 2 * lambda;
	double K2 = (PI / 2) * local_Cs * lambda * lambda * lambda;
	double K3 = atan(Q0/sqrt(1-Q0*Q0));

	_offsetCtf(local_Cs, lambda, rad_azimuth, defocus_average, defocus_deviation,
			   K1, K2, K3, mdt, particle);
}


OriginalBasis::OriginalBasis()
	:   AberrationBasis(5)
{}

void OriginalBasis::getBasisValues(double x, double y, double *dest)
{
	dest[0] = 1.0; // phase shift
	dest[1] = x*x + y*y; // defocus
	dest[2] = x*x - y*y; // oblique astigmatism
	dest[3] = x*y; // vertical astigmatism
	dest[4] = (x*x + y*y)*(x*x + y*y); // primary spherical
}

void OriginalBasis::_offsetCtf(
		double local_Cs, double lambda,
		double rad_azimuth, double defocus_average, double defocus_deviation,
		double K1, double K2, double K3, MetaDataTable &mdt, int particle)
{
	/* from ctf.h:

		   RFLOAT u2 = X * X + Y * Y;
		   RFLOAT u4 = u2 * u2;
		   RFLOAT deltaf = defocus_average + defocus_deviation*cos(2*(atan2(Y, X) - rad_azimuth))

		   argument = K1 * deltaf * u2 + K2 * u4 - K5 - K3

				K1 = PI / 2 * 2 * lambda;
				K2 = PI / 2 * local_Cs * lambda * lambda * lambda;
				K3 = atan(Q0/sqrt(1-Q0*Q0));
				K5 = DEG2RAD(phase_shift);

				local_Cs = Cs * 1e7;

	   astigmatism/defocus:

		   K1 * deltaf * u2
		   = K1 * defocus_average * u2 + defocus_deviation * K1 * cos(2*(phi - rad_azimuth)) * u2
		   = K1 * defocus_average * u2 + defocus_deviation * K1 * cos(2*phi - 2*rad_azimuth) * u2
		   = K1 * defocus_average * u2 + defocus_deviation * K1 * [cos(2*phi) cos(2*rad_azimuth) + sin(2*phi) sin(2*rad_azimuth)] * u2
		   = K1 * defocus_average * u2 + defocus_deviation * K1 * [(cos²(phi) - sin²(phi)) cos(2*rad_azimuth) + 2 sin(phi) cos(phi) sin(2*rad_azimuth)] * u2
		   = K1 * defocus_average * u2 + defocus_deviation * K1 * [(X² - Y²) cos(2*rad_azimuth) + 2 Y X sin(2*rad_azimuth)]
		   = b1 (X² + Y²) + b2 (X² - Y²) + b3 (XY)

		   where:  b1 = K1 * defocus_average
				   b2 = K1 * defocus_deviation * cos(2*rad_azimuth)
				   b3 = 2 * K1 * defocus_deviation * sin(2*rad_azimuth)

				   <=>

				   defocus_average = b1 / (PI lambda)
				   defocus_deviation = sqrt(b2² + b3²/4)/(PI lambda)
				   rad_azimuth = atan2(b3/2, b2) / 2                        */

	double b1 = K1 * defocus_average + coefficients[1];
	double b2 = K1 * defocus_deviation * cos(2*rad_azimuth) + coefficients[2];
	double b3 = 2 * K1 * defocus_deviation * sin(2*rad_azimuth) + coefficients[3];

	double new_defocus_average = b1 / (PI * lambda);
	double new_defocus_deviation = sqrt(b2*b2 + b3*b3/4)/(PI*lambda);
	double new_rad_azimuth = atan2(b3/2.0, b2) / 2.0;

	double azimuthal_angle = RAD2DEG(new_rad_azimuth);
	double DeltafU = -new_defocus_average - new_defocus_deviation;
	double DeltafV = new_defocus_deviation - new_defocus_average;

	/*     spherical aberration:

		   K2 * u4 = b4 * u4
			<=>
		   PI / 2 * local_Cs * lambda * lambda * lambda = b4;
			<=>
		   local_Cs = 2 * b4 / (PI lambda³)
			<=>
		   Cs = 1e-7 * 2 * b4 / (PI lambda³)            */

	double b4 = PI * lambda*lambda*lambda * local_Cs / 2.0 + coefficients[4];
	double Cs = 1e-7 * 2.0 * b4 / (PI * lambda*lambda*lambda);

	/*     phase shift / amp. contrast:

			   K3 = atan(Q0/sqrt(1-Q0*Q0))
				<=>
			   Q0/sqrt(1-Q0²) = tan(K3)
				<=>
			   Q0² = (1 - Q0²) * tan²(K3)
				<=>
			   Q0²(1 + tan²K3) = tan²K3
				<=>
			   Q0 = sqrt(tan²K3/(1 + tan²K3))
	*/

	double b0 = K3 - coefficients[0];

	if (b0 < 0)
	{
		double phase_shift;
		if (!mdt.getValue(EMDL_CTF_PHASESHIFT, phase_shift, particle)) phase_shift = 0;

		phase_shift = phase_shift - RAD2DEG(coefficients[0]);

		mdt.setValue(EMDL_CTF_PHASESHIFT, phase_shift, particle);
	}
	else
	{
		double t0 = tan(b0);
		double Q0 = sqrt(t0*t0/(1 + t0*t0));

		mdt.setValue(EMDL_CTF_Q0, Q0, particle);
	}

	mdt.setValue(EMDL_CTF_DEFOCUSU, DeltafU, particle);
	mdt.setValue(EMDL_CTF_DEFOCUSV, DeltafV, particle);
	mdt.setValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, particle);
	mdt.setValue(EMDL_CTF_CS, Cs, particle);

	//std::cout << DeltafU << ", " << DeltafV << " @ " << azimuthal_angle << "°, " << Cs << ", " << Q0 << "\n\n";

}

EvenData &EvenData::operator+=(const EvenData& d)
{
	Axx += d.Axx;
	Axy += d.Axy;
	Ayy += d.Ayy;
	bx += d.bx;
	by += d.by;

	return *this;
}

EvenData &EvenData::operator*=(double d)
{
	Axx *= d;
	Axy *= d;
	Ayy *= d;
	bx *= d;
	by *= d;

	return *this;
}

void EvenData::write(const RawImage<EvenData> &data, std::string filename)
{
	const int sh = data.xdim;
	const int s  = data.ydim;
	const int fc = data.zdim;

	BufferedImage<float>
		Axx(sh,s,fc),
		Axy(sh,s,fc),
		Ayy(sh,s,fc),
		bx(sh,s,fc),
		by(sh,s,fc);

	for (int f = 0; f < fc; f++)
	for (int y = 0; y <  s; y++)
	for (int x = 0; x < sh; x++)
	{
		Axx(x,y,f) = data(x,y,f).Axx;
		Axy(x,y,f) = data(x,y,f).Axy;
		Ayy(x,y,f) = data(x,y,f).Ayy;
		bx(x,y,f)  = data(x,y,f).bx;
		by(x,y,f)  = data(x,y,f).by;
	}

	Axx.write(filename + "_Axx.mrc");
	Axy.write(filename + "_Axy.mrc");
	Ayy.write(filename + "_Ayy.mrc");
	bx.write( filename + "_bx.mrc");
	by.write( filename + "_by.mrc");
}

BufferedImage<EvenData> EvenData::read(std::string filename)
{
	BufferedImage<double>
		Axx,
		Axy,
		Ayy,
		bx,
		by;

	Axx.read(filename + "_Axx.mrc");
	Axy.read(filename + "_Axy.mrc");
	Ayy.read(filename + "_Ayy.mrc");
	bx.read( filename + "_bx.mrc");
	by.read( filename + "_by.mrc");

	const int sh = Axx.xdim;
	const int s  = Axx.ydim;
	const int fc = Axx.zdim;

	BufferedImage<EvenData> data(sh,s,fc);

	for (int f = 0; f < fc; f++)
	for (int y = 0; y <  s; y++)
	for (int x = 0; x < sh; x++)
	{
		data(x,y,f).Axx = Axx(x,y,f);
		data(x,y,f).Axy = Axy(x,y,f);
		data(x,y,f).Ayy = Ayy(x,y,f);
		data(x,y,f).bx  = bx(x,y,f);
		data(x,y,f).by  = by(x,y,f);
	}

	return data;
}

OddData &OddData::operator+=(const OddData& d)
{
	a += d.a;
	b += d.b;

	return *this;
}

OddData &OddData::operator*=(double d)
{
	a *= d;
	b *= d;

	return *this;
}

void OddData::write(const RawImage<OddData> &data, std::string filename)
{
	const int sh = data.xdim;
	const int s  = data.ydim;
	const int fc = data.zdim;

	BufferedImage<float>
		a(sh,s,fc),
		b_real(sh,s,fc),
		b_imag(sh,s,fc);

	for (int f = 0; f < fc; f++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		a(x,y,f) = data(x,y,f).a;
		b_real(x,y,f) = data(x,y,f).b.real;
		b_imag(x,y,f) = data(x,y,f).b.imag;
	}

	a.write(filename + "_a.mrc");
	b_real.write(filename + "_b_real.mrc");
	b_imag.write(filename + "_b_imag.mrc");
}

BufferedImage<OddData> OddData::read(std::string filename)
{
	BufferedImage<double>
		a,
		b_real,
		b_imag;

	a.read(filename + "_a.mrc");
	b_real.read(filename + "_b_real.mrc");
	b_imag.read(filename + "_b_imag.mrc");

	const int sh = a.xdim;
	const int s  = a.ydim;
	const int fc = a.zdim;

	BufferedImage<OddData> data(sh,s,fc);

	for (int f = 0; f < fc; f++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		data(x,y,f).a = a(x,y,f);
		data(x,y,f).b.real = b_real(x,y,f);
		data(x,y,f).b.imag = b_imag(x,y,f);
	}

	return data;
}
