#include "Zernike_helper.h"
#include <src/jaz/optimization/nelder_mead.h>
#include <src/ctf.h>
#include <src/jaz/util/zio.h>

using namespace gravis;


BufferedImage<double> ZernikeHelper::computeEvenZernike(
		int s, double pixelSize, const d2Matrix& mag, int n_max)
{
	const int cc = Zernike::numberOfEvenCoeffs(n_max);
	const double as = (double)s * pixelSize;
	const int sh = s/2 + 1;

	BufferedImage<double> basis(sh,s,cc);

	for (int c = 0; c < cc; c++)
	{
		int m, n;
		Zernike::evenIndexToMN(c, m, n);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			const double xx0 = x/as;
			const double yy0 = y < s/2? y/as : (y-s)/as;

			const double xx = mag(0,0) * xx0 + mag(0,1) * yy0;
			const double yy = mag(1,0) * xx0 + mag(1,1) * yy0;

			basis(x,y,c) = Zernike::Z_cart(m, n, xx, yy);
		}
	}

	return basis;
}

std::vector<double> ZernikeHelper::fitEvenZernike(
		const BufferedImage<double>& phase,
		const BufferedImage<double>& weight,
		double pixelSize, const d2Matrix& mag, int n_max,
		BufferedImage<double>* fit)
{
	const int w = phase.xdim;
	const int h = phase.ydim;
	const int cc = Zernike::numberOfEvenCoeffs(n_max);

	BufferedImage<double> basis = computeEvenZernike(h, pixelSize, mag, n_max);

	std::vector<double> out = fitBasisLin(phase, weight, basis);

	if (fit != 0)
	{
		*fit = BufferedImage<double>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			(*fit)(x,y) = 0.0;
			
			for (int c = 0; c < cc; c++)
			{
				(*fit)(x,y) += out[c] * basis(x,y,c);
			}
		}
	}

	return out;
}

std::vector<double> ZernikeHelper::optimiseEvenZernike(
	const BufferedImage<dComplex>& xy,
	const BufferedImage<Tensor2x2<double>>& A,
	double pixelSize, const d2Matrix& mag, int n_max,
	const std::vector<double>& coeffs,
	BufferedImage<double>* fit)
{
	const int w = xy.xdim;
	const int h = xy.ydim;
	const int cc = Zernike::numberOfEvenCoeffs(n_max);

	BufferedImage<double> basis = computeEvenZernike(h, pixelSize, mag, n_max);

	std::vector<double> opt = optimiseBasis(xy, A, basis, coeffs);

	if (fit != 0)
	{
		*fit = BufferedImage<double>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(x,y) += opt[c] * basis(x,y,c);
			}
		}
	}

	return opt;
}



BufferedImage<double> ZernikeHelper::computeOddZernike(
		int s, double pixelSize, const d2Matrix& mag, int n_max)
{
	const int cc = Zernike::numberOfOddCoeffs(n_max);
	const double as = (double)s * pixelSize;
	const int sh = s/2 + 1;
	
	BufferedImage<double> basis(sh,s,cc);

	for (int c = 0; c < cc; c++)
	{
		int m, n;

		Zernike::oddIndexToMN(c, m, n);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			const double xx0 = x/as;
			const double yy0 = y < s/2? y/as : (y-s)/as;

			const double xx = mag(0,0) * xx0 + mag(0,1) * yy0;
			const double yy = mag(1,0) * xx0 + mag(1,1) * yy0;

			basis(x,y,c) = Zernike::Z_cart(m, n, xx, yy);
		}
	}

	return basis;
}

std::vector<double> ZernikeHelper::fitOddZernike(
	const BufferedImage<dComplex>& xy,
	const BufferedImage<double>& weight,
	double pixelSize, const d2Matrix& mag,
	int n_max,
	BufferedImage<double>* fit)
{
	const int w = xy.xdim;
	const int h = xy.ydim;
	const int cc = Zernike::numberOfOddCoeffs(n_max);

	BufferedImage<double> basis = ZernikeHelper::computeOddZernike(h, pixelSize, mag, n_max);

	std::vector<double> out = fitBasisLin(xy, weight, basis);

	if (fit != 0)
	{
		*fit = BufferedImage<double>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(x,y) += out[c] * basis(x,y,c);
			}
		}
	}

	return out;
}

std::vector<double> ZernikeHelper::optimiseOddZernike(
	const BufferedImage<dComplex> &xy,
	const BufferedImage<double> &weight,
	double pixelSize, const d2Matrix& mag, int n_max,
	const std::vector<double> &coeffs,
	BufferedImage<double> *fit)
{
	const int w = xy.xdim;
	const int h = xy.ydim;
	const int cc = Zernike::numberOfOddCoeffs(n_max);

	BufferedImage<double> basis = computeOddZernike(h, pixelSize, mag, n_max);

	std::vector<double> opt = optimiseBasis(xy, weight, basis, coeffs);

	if (fit != 0)
	{
		*fit = BufferedImage<double>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(x,y) += opt[c] * basis(x,y,c);
			}
		}
	}

	return opt;
}



std::vector<double> ZernikeHelper::fitBasisLin(
	const BufferedImage<dComplex>& xy,
	const BufferedImage<double>& weight,
	const BufferedImage<double>& basis)
{
	const int w = xy.xdim;
	const int h = xy.ydim;

	BufferedImage<double> phase(w,h);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		phase(x,y) = xy(x,y).arg();
	}

	return fitBasisLin(phase, weight, basis);
}

std::vector<double> ZernikeHelper::fitBasisLin(
	const BufferedImage<double>& phase,
	const BufferedImage<double>& weight,
	const BufferedImage<double>& basis)
{
	const int cc = basis.zdim;
	const int w = phase.xdim;
	const int h = phase.ydim;

	Matrix2D<double> A(cc,cc);
	Matrix1D<double> b(cc);

	for (int c1 = 0; c1 < cc; c1++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			b(c1) += weight(x,y) * basis(x,y,c1) * phase(x,y);
		}

		for (int c2 = c1; c2 < cc; c2++)
		{
			for (int y = 0; y < h; y++)
			for (int x = 0; x < w; x++)
			{
				A(c1,c2) += weight(x,y) * basis(x,y,c1) * basis(x,y,c2);
			}
		}

		for (int c2 = 0; c2 < c1; c2++)
		{
			A(c1,c2) = A(c2,c1);
		}
	}

	const double tol = 1e-20;
	Matrix1D<RFLOAT> x(cc);
	solve(A, b, x, tol);

	std::vector<double> out(cc);

	for (int c = 0; c < cc; c++)
	{
		out[c] = x(c);
	}

	return out;
}

std::vector<double> ZernikeHelper::optimiseBasis(
		const BufferedImage<dComplex>& xy,
		const BufferedImage<double>& weight,
		const BufferedImage<double>& basis,
		const std::vector<double>& initial)
{
	BasisOptimisation prob(xy, weight, basis, false);

	std::vector<double> opt = NelderMead::optimize(
		initial, prob, 0.01, 0.000001, 100000, 1.0, 2.0, 0.5, 0.5, false);

	return opt;
}

std::vector<double> ZernikeHelper::optimiseBasis(
		const BufferedImage<dComplex>& xy,
		const BufferedImage<Tensor2x2<double>>& A,
		const BufferedImage<double>& basis,
		const std::vector<double>& initial)
{
	AnisoBasisOptimisation prob(xy, A, basis, false);

	std::vector<double> opt = NelderMead::optimize(
		initial, prob, 0.01, 0.000001, 100000, 1.0, 2.0, 0.5, 0.5, false);

	return opt;
}

ZernikeHelper::OldCtfBasis ZernikeHelper::paramsToCtf(
		const std::vector<double> coeffs, double kV)
{
	OldCtfBasis out;
	
	const double local_kV = kV * 1e3;
	const double lambda = 12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));
	
	const int cc = coeffs.size();
	const double z0 = cc > 0? coeffs[0] : 0.0;
	const double z1 = cc > 1? coeffs[1] : 0.0;
	const double z2 = cc > 2? coeffs[2] : 0.0;
	const double z3 = cc > 3? coeffs[3] : 0.0;
	const double z6 = cc > 6? coeffs[6] : 0.0;
		
	const double K2_pix = 6.0 * z6;
	const double K3 = z2 - z0 - z6;

	std::cout << "cc = " << cc << std::endl;
	std::cout << "z0 = " << z0 << std::endl;
	std::cout << "z2 = " << z2 << std::endl;
	std::cout << "z6 = " << z6 << std::endl;
	std::cout << "K3 = " << K3 << std::endl;

	const double a = 2.0 * z2 - 6.0 * z6;
	const double b = z3;
	const double c = z1;
	
	double ev0, ev1, cs, sn;

	dsyev2(a+b, c, a-b, &ev0, &ev1, &cs, &sn);

	const double defU_pix = ev0;
	const double defV_pix = ev1;
	
	const double azimuth_deg = RAD2DEG(atan2(sn,cs));
	
	const double defConv = 1.0 / (PI * lambda);
	
	
	const double tanK3 = tan(K3);
	const double tan2K3 = tanK3*tanK3;
	
	out.defocusU = -defConv * defU_pix;
	out.defocusV = -defConv * defV_pix;
	out.astigAzimuth_deg = azimuth_deg;
	
	out.Cs = 2.0 * 1e-7 * K2_pix / (PI * lambda*lambda*lambda);
		
	if (K3 > 0.0 && K3 < atan(sqrt(0.5)))
	{
		out.Q0 = K3 > 0.0? sqrt(tan2K3 / (1.0 + tan2K3)) : -sqrt(tan2K3 / (1.0 + tan2K3));
		out.phaseShift_deg = 0.0;
	}
	else
	{
		out.Q0 = 0.0;
		out.phaseShift_deg = RAD2DEG(K3);
	}
	
	return out;
}

std::vector<double> ZernikeHelper::ctfToParams(const CTF &ctf)
{
	const double local_kV = ctf.kV * 1e3;
	const double rad_azimuth = DEG2RAD(ctf.azimuthal_angle);

	const double lambda = 12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));

	const double K1 = (PI / 2) * 2 * lambda;
	const double K3 = atan(ctf.Q0 / sqrt(1.0 - ctf.Q0 * ctf.Q0));
	const double K5 = DEG2RAD(ctf.phase_shift);

	const double sin_az = sin(rad_azimuth);
	const double cos_az = cos(rad_azimuth);

	d2Matrix Q(cos_az, sin_az, -sin_az, cos_az);
	d2Matrix Qt(cos_az, -sin_az, sin_az, cos_az);
	d2Matrix D(-ctf.DeltafU, 0.0, 0.0, -ctf.DeltafV);

	d2Matrix A = K1 * Qt * D * Q;

	const double Axx = A(0,0);
	const double Axy = A(0,1);
	const double Ayy = A(1,1);
	
	std::vector<double> out(8);
		
	const double a = (Axx + Ayy) / 2.0;
	const double b = (Axx - Ayy) / 2.0;
	const double c = Axy;
	
	const double K2_pix = ctf.Cs * (PI * lambda*lambda*lambda) / (2.0 * 1e-7);
	
	out[0] = -K5 - K3 + (a + K2_pix) / 2.0 - K2_pix / 6.0;
	out[1] = c;
	out[2] = (a + K2_pix) / 2.0;
	out[3] = b;
	out[6] = K2_pix / 6.0;
	
	return out;
}

void ZernikeHelper::testConversion()
{
	const int iters = 10;
	const double pixelSize = 1.3;
	const double kV = 300.0;
	const int s = 512;
	const int sh = s/2 + 1;
	
	const double local_kV = kV * 1e3;
	const double lambda = 12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));
	
	
	BufferedImage<double> basis = computeEvenZernike(s, pixelSize, d2Matrix(), 4);
	basis.write("debug/basis.mrc");
	
	
	for (int it = 0; it < iters; it++)
	{
		std::vector<double> coeffs(8, 0.0);
		
		coeffs[0] = -0.5 * (rand() / (double) RAND_MAX - 0.5);
		coeffs[1] = rand() / (double) RAND_MAX;
		coeffs[2] = 3.0 * rand() / (double) RAND_MAX;
		coeffs[3] = rand() / (double) RAND_MAX;
		coeffs[6] = 0.3 * rand() / (double) RAND_MAX;
		
		std::cout 
			<< coeffs[0] << " \t" 
			<< coeffs[1] << " \t" 
			<< coeffs[2] << " \t" 
			<< coeffs[3] << " \t" 
			<< coeffs[6] << std::endl; 
		
		
		BufferedImage<double> comb0(sh,s);
		comb0.fill(0.0);
		
		for (int i = 0; i < 8; i++)
		{
			comb0 += coeffs[i] * basis.getSliceRef(i);
		}
		
		comb0.write("debug/"+ZIO::itoa(it)+"_Zernike.mrc");
		
		//std::cout << "comb0(0,0) = " << comb0(0,0) << std::endl;
		
		
		
		OldCtfBasis ob = paramsToCtf(coeffs, kV);
				
		CTF ctf;
		ctf.setValues(ob.defocusU, ob.defocusV, ob.astigAzimuth_deg, kV, ob.Cs, ob.Q0,
					  0.0, 1.0, ob.phaseShift_deg);
						
		BufferedImage<double> comb1(sh,s);
		ctf.drawGamma(s, s, pixelSize, &comb1[0]);
		
		comb1.write("debug/"+ZIO::itoa(it)+"_CTF.mrc");
		
		double L2 = 0.0;
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double d = comb1(x,y) - comb0(x,y);
			L2 += d*d;
		}
		
		//std::cout << "L2: " << L2 << std::endl;
		
		(comb1 - comb0).write("debug/"+ZIO::itoa(it)+"_diff.mrc");
		
		std::vector<double> coeffs2 = ctfToParams(ctf);
		
		std::cout 
			<< coeffs2[0] << " \t" 
			<< coeffs2[1] << " \t" 
			<< coeffs2[2] << " \t" 
			<< coeffs2[3] << " \t" 
			<< coeffs2[6] << std::endl; 
		
		double L2c = 0.0;
		
		for (int i = 0; i < 8; i++)
		{
			if (i > 3 && i != 6) continue;			
			const double d = coeffs2[i] - coeffs[i];
			L2c += d*d;
		}
		
		//std::cout << "L2 coeff: " << L2c << std::endl;
		
	}
}



ZernikeHelper::BasisOptimisation::BasisOptimisation(
		const BufferedImage<dComplex> &xy,
		const BufferedImage<double> &weight,
		const BufferedImage<double> &basis,
		bool L1)
:	w(xy.xdim),
	h(xy.ydim),
	cc(basis.zdim),
	xy(xy),
	weight(weight),
	basis(basis),
	L1(L1)
{
}

double ZernikeHelper::BasisOptimisation::f(const std::vector<double>& x, void* tempStorage) const
{
	BufferedImage<double>& recomb = *((BufferedImage<double>*)tempStorage);
	recomb.fill(0.0);

	for (int c  = 0; c < cc; c++)
	for (int yp = 0; yp < h; yp++)
	for (int xp = 0; xp < w; xp++)
	{
		recomb(xp,yp) += x[c] * basis(xp,yp,c);
	}

	double sum = 0.0;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		Complex zPred(cos(recomb(x,y)), sin(recomb(x,y)));
		sum += weight(x,y) * (zPred - xy(x,y)).norm();
	}

	return sum;
}

void* ZernikeHelper::BasisOptimisation::allocateTempStorage() const
{
	return new BufferedImage<double>(w,h);
}

void ZernikeHelper::BasisOptimisation::deallocateTempStorage(void* ts)
{
	delete (BufferedImage<double>*)ts;
}






ZernikeHelper::AnisoBasisOptimisation::AnisoBasisOptimisation(
		const BufferedImage<dComplex>& xy,
		const BufferedImage<Tensor2x2<double>>& A,
		const BufferedImage<double>& basis,
		bool L1)
:	w(xy.xdim),
	h(xy.ydim),
	cc(basis.zdim),
	xy(xy),
	A(A),
	basis(basis),
	L1(L1)
{
}

double ZernikeHelper::AnisoBasisOptimisation::f(const std::vector<double>& x, void* tempStorage) const
{
	BufferedImage<double>& recomb = *((BufferedImage<double>*)tempStorage);
	recomb.fill(0.0);

	for (int c  = 0; c < cc; c++)
	for (int yp = 0; yp < h; yp++)
	for (int xp = 0; xp < w; xp++)
	{
		recomb(xp,yp) += x[c] * basis(xp,yp,c);
	}

	double sum = 0.0;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		d2Vector e(cos(recomb(x,y)) - xy(x,y).real, sin(recomb(x,y)) - xy(x,y).imag);

		sum +=
					A(x,y).xx * e.x * e.x
			+ 2.0 * A(x,y).xy * e.x * e.y
			+       A(x,y).yy * e.y * e.y;
	}

	return sum;
}

void* ZernikeHelper::AnisoBasisOptimisation::allocateTempStorage() const
{
	return new BufferedImage<double>(w,h);
}

void ZernikeHelper::AnisoBasisOptimisation::deallocateTempStorage(void *ts)
{
	delete (BufferedImage<double>*)ts;
}





ZernikeHelper::MultiAnisoBasisOptimisation::MultiAnisoBasisOptimisation(const BufferedImage<dComplex>& xy,
		const BufferedImage<Tensor2x2<double>>& A,
		const BufferedImage<double>& basis,
		double lambda,
		bool L1)
:
	w(xy.xdim),
	h(xy.ydim),
	fc(xy.zdim),
	cc(basis.zdim),
	lambda(lambda),
	xy(xy),
	A(A),
	basis(basis),
	L1(L1)
{
}

double ZernikeHelper::MultiAnisoBasisOptimisation::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
{
	BufferedImage<double> recomb(w,h);
	recomb.fill(0.0);

	for (int i = 0; i < gradDest.size(); i++)
	{
		gradDest[i] = 0.0;
	}

	double sum = 0.0;

	for (int f = 0; f < fc; f++)
	{
		recomb.fill(0.0);

		for (int c  = 0; c < cc; c++)
		for (int yy = 0; yy < h; yy++)
		for (int xx = 0; xx < w; xx++)
		{
			recomb(xx,yy) += x[3*f + c + 3] * basis(xx,yy,c);
		}

		for (int yy = 0; yy < h; yy++)
		for (int xx = 0; xx < w; xx++)
		{
			const dComplex rot(cos(recomb(xx,yy)), sin(recomb(xx,yy)));
			const dComplex e(rot.real - xy(xx,yy,f).real, rot.imag - xy(xx,yy,f).imag);

			sum +=
						A(xx,yy,f).xx * e.real * e.real
				+ 2.0 * A(xx,yy,f).xy * e.real * e.imag
				+       A(xx,yy,f).yy * e.imag * e.imag;

			// drecomb_dx[3*f + c + 3] = basis(xp,yp,c)

			const dComplex drot_drecomb(-sin(recomb(xx,yy)), cos(recomb(xx,yy)));

			// de_drot = 1, 1

			const dComplex dsum_de(
				2.0 * (A(xx,yy,f).xx * e.real + A(xx,yy,f).xy * e.imag),
				2.0 * (A(xx,yy,f).yy * e.imag + A(xx,yy,f).xy * e.real));

			const double dsum_drecomb =
					  dsum_de.real * drot_drecomb.real
					+ dsum_de.imag * drot_drecomb.imag;

			for (int c  = 0; c < cc; c++)
			{
				gradDest[3*f + c + 3] += dsum_drecomb * basis(xx,yy,c);
			}
		}
	}

	for (int i = 3; i < gradDest.size(); i++)
	{
		const double d = x[i] - x[i%3];
		const double q = lambda * w * h / 1e6;

		sum += q * d * d;

		gradDest[i]   += 2.0 * q * d;
		gradDest[i%3] -= 2.0 * q * d;
	}

	return sum;
}

double ZernikeHelper::MultiAnisoBasisOptimisation::f(const std::vector<double>& x, void* tempStorage) const
{
	BufferedImage<double>& recomb = *((BufferedImage<double>*)tempStorage);

	double sum = 0.0;

	for (int f = 0; f < fc; f++)
	{
		recomb.fill(0.0);

		for (int c  = 0; c < cc; c++)
		for (int yp = 0; yp < h; yp++)
		for (int xp = 0; xp < w; xp++)
		{
			recomb(xp,yp) += x[3*f + c + 3] * basis(xp,yp,c);
		}

		for (int yy = 0; yy < h; yy++)
		for (int xx = 0; xx < w; xx++)
		{
			d2Vector e(cos(recomb(xx,yy)) - xy(xx,yy,f).real, sin(recomb(xx,yy)) - xy(xx,yy,f).imag);

			sum +=
						A(xx,yy,f).xx * e.x * e.x
				+ 2.0 * A(xx,yy,f).xy * e.x * e.y
				+       A(xx,yy,f).yy * e.y * e.y;
		}
	}

	return sum;
}

void* ZernikeHelper::MultiAnisoBasisOptimisation::allocateTempStorage() const
{
	return new BufferedImage<double>(w,h);
}

void ZernikeHelper::MultiAnisoBasisOptimisation::deallocateTempStorage(void *ts)
{
	delete (BufferedImage<double>*)ts;
}
