#include "modular_ctf_optimisation.h"
#include <src/jaz/ctf/magnification_helper.h>
#include <src/ctf.h>

using namespace gravis;



ModularCtfOptimisation::ModularCtfOptimisation(
	MetaDataTable &mdt, 
	ObservationModel* obsModel, 
	const std::vector<Image<Complex>>& obs, 
	const std::vector<Image<Complex>>& pred, 
	const std::vector<Image<RFLOAT>>& frqWghByGroup,
	std::string modeStr)
	:
	mdt(mdt),
	obsModel(obsModel),
	obs(obs),
	pred(pred),
	particle_count(mdt.numberOfObjects()),
	frqWghByGroup(frqWghByGroup)
{
	initialValues.resize(CtfParamCount * particle_count);
	
	for (int p = 0; p < particle_count; p++)
	{
		CTF ctf;
		ctf.readByGroup(mdt, obsModel, p);
		
		std::vector<double> K = ctf.getK();
		
		const double Axx = K[1] * ctf.getAxx();
		const double Axy = K[1] * ctf.getAxy();
		const double Ayy = K[1] * ctf.getAyy();
				
		const double dz = (Axx + Ayy) / 2.0;
		const double bfac = ctf.Bfac;
		const double kfac = ctf.scale;
				
		initialValues[CtfParamCount * p + Phase]               = -K[5] - K[3];
		initialValues[CtfParamCount * p + Defocus]             = dz;
		initialValues[CtfParamCount * p + Astigmatism1]        = Axx - dz;
		initialValues[CtfParamCount * p + Astigmatism2]        = Axy;
		initialValues[CtfParamCount * p + SphericalAberration] = K[2];
		initialValues[CtfParamCount * p + BFactor]             = bfac;
		initialValues[CtfParamCount * p + ScaleFactor]         = kfac;
		
		/*if (p == 0)
		{
			std::cout << "p0_1: " << ctf.DeltafU << " " << ctf.DeltafV << " " << ctf.azimuthal_angle << "\n";
		}*/
	}
	
	modes = decodeModes(modeStr);
			
	int currOff = 0;
	
	for (int i = 0; i < CtfParamCount; i++)
	{
		switch (modes[i])
		{
			case PerParticle:
			{
				paramOffset[i] = currOff;
				paramParticleStep[i] = 1;
				
				currOff += particle_count;
				
				break;
			}
				
			case PerMicrograph:
			{
				paramOffset[i] = currOff;
				paramParticleStep[i] = 0;
				
				currOff += 1;
				
				break;
			}
				
			case Fixed:
			{
				paramOffset[i] = 0;
				paramParticleStep[i] = 0;
				
				break;
			}
		}
	}
	
	param_count = currOff;
	
	angpix = obsModel->getPixelSizes();
	
	if (obsModel->hasEvenZernike)
	{
		std::vector<int> ogPres = obsModel->getOptGroupsPresent_zeroBased(mdt);
		
		aberrationByGroup.resize(obsModel->numberOfOpticsGroups());
		
		for (int ogpi = 0; ogpi < ogPres.size(); ogpi++)
		{
			const int og = ogPres[ogpi];			
			aberrationByGroup[og] = obsModel->getGammaOffset(og, obsModel->getBoxSize(og));
		}
	}
}

double ModularCtfOptimisation::f(const std::vector<double> &x, void *tempStorage) const
{
	double out = 0.0;
	
	for (int p = 0; p < particle_count; p++)
	{
		const double ph = readParam(Phase, x, p);
		const double dz = readParam(Defocus, x, p);
		const double a1 = readParam(Astigmatism1, x, p);
		const double a2 = readParam(Astigmatism2, x, p);
		const double cs = readParam(SphericalAberration, x, p);
		const double bf = readParam(BFactor, x, p);
		const double kf = readParam(ScaleFactor, x, p);
		
		const int s = obs[p].data.ydim;
		const int sh = s/2 + 1;
		
		const int og = obsModel->getOpticsGroup(mdt, p);
		const double as = s * angpix[og];
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double wght = (x == 0)? 1.0 : 2.0;
			
			const double xx = x/as;
			const double yy = (y <= s/2)? y/as : (y - s)/as;
			
			double xu, yu;
			
			if (obsModel->hasMagMatrices)
			{
				const Matrix2D<RFLOAT>& M = obsModel->getMagMatrix(og);
				
				xu = M(0,0) * xx + M(0,1) * yy;
				yu = M(1,0) * xx + M(1,1) * yy;
			}
			else
			{
				xu = xx;
				yu = yy;
			}
			
			const double u2 = xu * xu + yu * yu;
			const double u4 = u2 * u2;
			
			const double gammaOffset = obsModel->hasEvenZernike? aberrationByGroup[og](y,x) : 0.0;
			const double freqWgh = frqWghByGroup[og](y,x);
	
			const double gamma = ph 
					+ (dz + a1) * xu * xu + 2.0 * a2 * xu * yu + (dz - a1) * yu * yu 
					+ cs * u4 + gammaOffset;
			
			const double ctfVal = -kf * exp(-bf * u2 / 4.0) * sin(gamma);
			
			out += freqWgh * wght * ( obs[p](y,x) - ctfVal * pred[p](y,x) ).norm();
		}
	}
	
	return out;
}

void ModularCtfOptimisation::grad(
		const std::vector<double> &x, 
		std::vector<double> &gradDest, 
		void *tempStorage) const
{
	for (int i = 0; i < gradDest.size(); i++)
	{
		gradDest[i] = 0.0;
	}
	
	for (int p = 0; p < particle_count; p++)
	{
		const double ph = readParam(Phase, x, p);
		const double dz = readParam(Defocus, x, p);
		const double a1 = readParam(Astigmatism1, x, p);
		const double a2 = readParam(Astigmatism2, x, p);
		const double cs = readParam(SphericalAberration, x, p);
		const double bf = readParam(BFactor, x, p);
		const double kf = readParam(ScaleFactor, x, p);
		
		const int s = obs[p].data.ydim;
		const int sh = s/2 + 1;
		
		const int og = obsModel->getOpticsGroup(mdt, p);
		const double as = s * angpix[og];
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double wght = (x == 0)? 1.0 : 2.0;
			
			const double xx = x/as;
			const double yy = (y <= s/2)? y/as : (y - s)/as;
			
			double xu, yu;
			
			if (obsModel->hasMagMatrices)
			{
				const Matrix2D<RFLOAT>& M = obsModel->getMagMatrix(og);
				
				xu = M(0,0) * xx + M(0,1) * yy;
				yu = M(1,0) * xx + M(1,1) * yy;
			}
			else
			{
				xu = xx;
				yu = yy;
			}
			
			const double u2 = xu * xu + yu * yu;
			const double u4 = u2 * u2;
			
			const double gammaOffset = obsModel->hasEvenZernike? aberrationByGroup[og](y,x) : 0.0;
			const double freqWgh = frqWghByGroup[og](y,x);
	
			const double gamma = ph 
					+ (dz + a1) * xu * xu + 2.0 * a2 * xu * yu + (dz - a1) * yu * yu 
					+ cs * u4 + gammaOffset;
			
			const double env = exp(-bf * u2 / 4.0);
			const double wave = -sin(gamma);
			const double ctfVal = kf * env * wave;
			
			const Complex z_obs = obs[p](y,x);
			const Complex z_pred = pred[p](y,x);
			const Complex z_err = ctfVal * z_pred - z_obs;
			
			const double dE_dCtfVal = 2.0 * freqWgh * wght * (
						z_err.real * z_pred.real + z_err.imag * z_pred.imag);
			
			const double dE_dGamma = dE_dCtfVal * kf * env * (-cos(gamma));
			
			const double dE_dPh = dE_dGamma;
			const double dE_dDz = dE_dGamma * (xu * xu + yu * yu);
			const double dE_dA1 = dE_dGamma * (xu * xu - yu * yu);
			const double dE_dA2 = dE_dGamma * 2.0 * a2 * xu * yu;
			const double dE_dCs = dE_dGamma * u4;
			const double dE_dBf = dE_dCtfVal * kf * wave * exp(-bf * u2 / 4.0) * (-u2/4.0);	
			const double dE_dKf = dE_dCtfVal * env * wave;
			
			if (modes[Phase] != Fixed)
			{
				gradDest[paramOffset[Phase] + p * paramParticleStep[Phase]] += dE_dPh;
			}
			
			if (modes[Defocus] != Fixed)
			{
				gradDest[paramOffset[Defocus] + p * paramParticleStep[Defocus]] += dE_dDz;
			}
			
			if (modes[Astigmatism1] != Fixed)
			{
				gradDest[paramOffset[Astigmatism1] + p * paramParticleStep[Astigmatism1]] += dE_dA1;
			}
			
			if (modes[Astigmatism2] != Fixed)
			{
				gradDest[paramOffset[Astigmatism2] + p * paramParticleStep[Astigmatism2]] += dE_dA2;
			}
			
			if (modes[SphericalAberration] != Fixed)
			{
				gradDest[paramOffset[SphericalAberration] + p * paramParticleStep[SphericalAberration]] += dE_dCs;
			}
			
			if (modes[BFactor] != Fixed)
			{
				gradDest[paramOffset[BFactor] + p * paramParticleStep[BFactor]] += dE_dBf;
			}
			
			if (modes[ScaleFactor] != Fixed)
			{
				gradDest[paramOffset[ScaleFactor] + p * paramParticleStep[ScaleFactor]] += dE_dKf;
			}
		}
	}
}

std::vector<double> ModularCtfOptimisation::encodeInitial()
{
	std::vector<double> out(param_count, 0.0);
	
	for (int cp = 0; cp < CtfParamCount; cp++)
	{
		switch (modes[cp])
		{
			case PerParticle:
			{
				for (int p = 0; p < particle_count; p++)
				{
					out[paramOffset[cp] + p] = initialValues[CtfParamCount * p + cp];
				}
				
				break;
			}
				
			case PerMicrograph:
			{
				out[paramOffset[cp]] = initialValues[cp]; // returning value of first particle
				
				break;
			}
				
			case Fixed:
			{				
				break;
			}
		}
	}
	
	return out;
}

void ModularCtfOptimisation::writeToTable(const std::vector<double> &x)
{
	for (int p = 0; p < particle_count; p++)
	{
		const double ph = readParam(Phase, x, p);
		const double dz = readParam(Defocus, x, p);
		const double a1 = readParam(Astigmatism1, x, p);
		const double a2 = readParam(Astigmatism2, x, p);
		const double cs = readParam(SphericalAberration, x, p);
		const double bf = readParam(BFactor, x, p);
		const double kf = readParam(ScaleFactor, x, p);
				
		CTF ctf;
		ctf.readByGroup(mdt, obsModel, p);
		
		std::vector<double> K = ctf.getK();		
		
		d2Matrix A(-(dz+a1)/ K[1],  -a2 / K[1],
				       -a2 / K[1],  -(dz-a1) / K[1]);
		
		RFLOAT defocusU, defocusV, angleDeg;		
		MagnificationHelper::matrixToPolar(A, defocusU, defocusV, angleDeg);
		
		RFLOAT local_kV = ctf.kV * 1e3;		
		const double lambda = 12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));
				
		/*if (p == 0)
		{
			std::cout << "p0_1: " << defocusU << " " << defocusV << " " << angleDeg << "\n";
		}*/
		
		// ph = -K[5] - K[3] = -DEG2RAD(phase_shift) - atan(Q0/sqrt(1-Q0*Q0));
		// => phase_shift = RAD2DEG(-ph - K[3])
		
		mdt.setValue(EMDL_CTF_PHASESHIFT, RAD2DEG(-ph - K[3]), p);	
		
		mdt.setValue(EMDL_CTF_DEFOCUSU, defocusU, p);
		mdt.setValue(EMDL_CTF_DEFOCUSV, defocusV, p);
		mdt.setValue(EMDL_CTF_DEFOCUS_ANGLE, angleDeg, p);
		
		// cs = K[2] = (PI / 2) * C_s * 1e7 * lambda^3
		// => C_s = cs * (2/PI) / (1e7 / lambda^3)
		
		mdt.setValue(EMDL_CTF_CS, 2 * cs / (1e7 * PI * lambda*lambda*lambda), p);
		
		mdt.setValue(EMDL_CTF_BFACTOR, bf, p);
		mdt.setValue(EMDL_CTF_SCALEFACTOR, kf, p);
	}
}

bool ModularCtfOptimisation::validateModeString(std::string mode)
{
	if (mode.length() != 5) return false;
	
	for (int i = 0; i < 5; i++)
	{
		if (mode[i] != 'p' && mode[i] != 'm' && mode[i] != 'f')
		{
			return false;
		}
	}
	
	return true;
}

std::vector<ModularCtfOptimisation::Mode> ModularCtfOptimisation::decodeModes(std::string s)
{
	std::vector<Mode> out(CtfParamCount);
	
	std::vector<Mode> charToMode(256, Fixed);
	
	charToMode[(unsigned char)'p'] = PerParticle;
	charToMode[(unsigned char)'m'] = PerMicrograph;
	
	out[Phase]                = charToMode[(unsigned char)s[0]];
	out[Defocus]              = charToMode[(unsigned char)s[1]];
	out[Astigmatism1]         = charToMode[(unsigned char)s[2]];
	out[Astigmatism2]         = charToMode[(unsigned char)s[2]];
	out[SphericalAberration]  = charToMode[(unsigned char)s[3]];
	out[BFactor]              = charToMode[(unsigned char)s[4]];
	out[ScaleFactor]          = charToMode[(unsigned char)s[4]];
	
	return out;
}
