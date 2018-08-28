#include "src/jaz/obs_model.h"
#include "src/jaz/stack_helper.h"
#include "src/jaz/filter_helper.h"
#include "src/jaz/Fourier_helper.h"
#include "src/jaz/ctf/tilt_helper.h"
#include "src/jaz/math/Zernike.h"
#include "src/jaz/vtk_helper.h"

#include <src/backprojector.h>

#include <set>
#include <omp.h>


void ObservationModel::loadSafely(
	std::string particlesFn, std::string opticsFn, 
	ObservationModel& obsModel, 
	MetaDataTable& particlesMdt, MetaDataTable& opticsMdt)
{
	particlesMdt.read(particlesFn);
	
	if (particlesMdt.getVersion() < 31000)
	{
		REPORT_ERROR_STR(particlesFn << " is from a pre-3.1 version of Relion. "
			<< "Please use relion_convert_star to generate an up-to-date file.");
	}

	if (!containsAllNeededColumns(particlesMdt))
	{
		REPORT_ERROR_STR(particlesFn << " does not contain all of the required columns ("
			<< "rlnOriginX, rlnOriginY, rlnAngleRot, rlnAngleTilt, rlnAnglePsi and rlnRandomSubset)");
	}
	
	opticsMdt.read(opticsFn);
		
	obsModel = ObservationModel(opticsMdt);	

	// read pixel sizes (and make sure they are all the same)	
	
	if (!obsModel.allPixelSizesIdentical())
	{
		REPORT_ERROR("ERROR: different pixel sizes detected. Please split your dataset by pixel size.");
	}
	
	// make sure all optics groups are defined
	
	std::vector<int> undefinedOptGroups = obsModel.findUndefinedOptGroups(particlesMdt);
	
	if (undefinedOptGroups.size() > 0)
	{
		std::stringstream sts;
		
		for (int i = 0; i < undefinedOptGroups.size(); i++)
		{
			sts << undefinedOptGroups[i];
			
			if (i < undefinedOptGroups.size()-1)
			{
				sts << ", ";
			}
		}
		
		REPORT_ERROR("ERROR: The following optics groups were not defined in "+
					 opticsFn+": "+sts.str());
	}
	
	// make sure the optics groups appear in the right order (and rename them if necessary)
	
	if (!obsModel.opticsGroupsSorted())
	{
		std::cerr << "   - Warning: the optics groups in " << opticsFn 
				  << " are not in the right order - renaming them now" << std::endl;
		
		obsModel.sortOpticsGroups(particlesMdt);
	}
}

ObservationModel::ObservationModel()
{
}

ObservationModel::ObservationModel(const MetaDataTable &opticsMdt)
:	opticsMdt(opticsMdt),
	angpix(opticsMdt.numberOfObjects()),
	lambda(opticsMdt.numberOfObjects()),
	Cs(opticsMdt.numberOfObjects())
{
	if (   !opticsMdt.containsLabel(EMDL_CTF_MAGNIFICATION)
	    || !opticsMdt.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE)
	    || !opticsMdt.containsLabel(EMDL_CTF_VOLTAGE)
	    || !opticsMdt.containsLabel(EMDL_CTF_CS))
	{
		REPORT_ERROR_STR("ERROR: not all necessary variables defined in _optics.star file: "
			<< "rlnMagnification, rlnDetectorPixelSize, rlnVoltage and rlnSphericalAberration.");
	}
	
	// symmetrical high-order aberrations:
	hasEvenZernike = opticsMdt.containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS);
	evenZernikeCoeffs = std::vector<std::vector<double> >(
			opticsMdt.numberOfObjects(), std::vector<double>(0));
	gammaOffset = std::vector<std::map<int,Image<RFLOAT> > >(opticsMdt.numberOfObjects());
	
	// antisymmetrical high-order aberrations:
	hasOddZernike = opticsMdt.containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS);		
	oddZernikeCoeffs = std::vector<std::vector<double> >(
			opticsMdt.numberOfObjects(), std::vector<double>(0));
	phaseCorr = std::vector<std::map<int,Image<Complex> > >(opticsMdt.numberOfObjects());
	
	const bool hasTilt = opticsMdt.containsLabel(EMDL_IMAGE_BEAMTILT_X)
	                  || opticsMdt.containsLabel(EMDL_IMAGE_BEAMTILT_Y);							 
	
	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		RFLOAT mag, dstep;
		opticsMdt.getValue(EMDL_CTF_MAGNIFICATION, mag, i);
		opticsMdt.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, i);
		angpix[i] = 10000 * dstep / mag;
		
		double kV;		
		opticsMdt.getValue(EMDL_CTF_VOLTAGE, kV, i);		
		double V = kV * 1e3;
		lambda[i] = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));
		
		opticsMdt.getValue(EMDL_CTF_CS, Cs[i], i);
		
		if (hasEvenZernike)
		{
			opticsMdt.getValue(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, evenZernikeCoeffs[i], i);
		}
		
		if (hasOddZernike)
		{
			opticsMdt.getValue(EMDL_IMAGE_ODD_ZERNIKE_COEFFS, oddZernikeCoeffs[i], i);
		}
		
		if (hasTilt)
		{
			double tx(0), ty(0);
			
			opticsMdt.getValue(EMDL_IMAGE_BEAMTILT_X, tx, i);
			opticsMdt.getValue(EMDL_IMAGE_BEAMTILT_Y, ty, i);			
			
			if (!hasOddZernike)
			{
				oddZernikeCoeffs[i] = std::vector<double>(6, 0.0);
			}
			
			TiltHelper::insertTilt(oddZernikeCoeffs[i], tx, ty, Cs[i], lambda[i]);
			
			hasOddZernike = true;
		}
	}
	
	// @TODO: make sure tilt is in opticsMDT, not in particlesMDT!
}

void ObservationModel::predictObservation(
        Projector& proj, const MetaDataTable& partMdt, long int particle,
		MultidimArray<Complex>& dest,
        bool applyCtf, bool shiftPhases, bool applyShift)
{
    const int s = proj.ori_size;
    const int sh = s/2 + 1;

    double xoff, yoff;

    partMdt.getValue(EMDL_ORIENT_ORIGIN_X, xoff, particle);
    partMdt.getValue(EMDL_ORIENT_ORIGIN_Y, yoff, particle);

    double rot, tilt, psi;

    Matrix2D<RFLOAT> A3D;
    partMdt.getValue(EMDL_ORIENT_ROT, rot, particle);
    partMdt.getValue(EMDL_ORIENT_TILT, tilt, particle);
    partMdt.getValue(EMDL_ORIENT_PSI, psi, particle);

    Euler_angles2matrix(rot, tilt, psi, A3D);

	if (dest.xdim != sh || dest.ydim != s)
	{
		dest.resize(s,sh);
	}
	
	dest.initZeros();

    proj.get2DFourierTransform(dest, A3D, false);
	
	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;
	
	if (applyShift)
	{
		shiftImageInFourierTransform(dest, dest, s, s/2 - xoff, s/2 - yoff);
	}

    if (applyCtf)
    {
        CTF ctf;
        ctf.readByGroup(partMdt, this, particle);
		
		Image<RFLOAT> ctfImg(sh,s);
		ctf.getFftwImage(ctfImg(), s, s, angpix[opticsGroup]);

		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			dest(y,x) *= ctfImg(y,x);
		}
    }

    if (shiftPhases && oddZernikeCoeffs.size() > opticsGroup 
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
    {
		const Image<Complex>& corr = getPhaseCorrection(opticsGroup, s);
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			dest(y,x) *= corr(y,x);
		}
    }
}

Image<Complex> ObservationModel::predictObservation(
        Projector& proj, const MetaDataTable& partMdt, long int particle,
        bool applyCtf, bool shiftPhases, bool applyShift)
{
    Image<Complex> pred;
	
	predictObservation(proj, partMdt, particle, pred.data, applyCtf, shiftPhases, applyShift);
    return pred;
}

std::vector<Image<Complex> > ObservationModel::predictObservations(
        Projector &proj, const MetaDataTable &partMdt, int threads,
        bool applyCtf, bool shiftPhases, bool applyShift)
{
    const int pc = partMdt.numberOfObjects();
    std::vector<Image<Complex> > out(pc);

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        out[p] = predictObservation(proj, partMdt, p, applyCtf, shiftPhases, applyShift);
    }

	return out;
}

void ObservationModel::demodulatePhase(
		const MetaDataTable& partMdt, long particle, MultidimArray<Complex>& obsImage)
{
	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;
	
	demodulatePhase(opticsGroup, obsImage);
}

void ObservationModel::demodulatePhase(
		int opticsGroup, MultidimArray<Complex>& obsImage)
{
	const int s = obsImage.ydim;
	const int sh = obsImage.xdim;
	
	if (oddZernikeCoeffs.size() > opticsGroup 
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
    {
		const Image<Complex>& corr = getPhaseCorrection(opticsGroup, s);
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			obsImage(y,x) *= corr(y,x).conj();
		}
    }
}

bool ObservationModel::allPixelSizesIdentical() const
{
	bool out = true;
	
	for (int i = 1; i < angpix.size(); i++)
	{
		if (angpix[i] != angpix[0])
		{
			out = false;
			break;
		}
	}
	
	return out;
}

double ObservationModel::angToPix(double a, int s, int opticsGroup) const
{
	return s * angpix[opticsGroup] / a;
}

double ObservationModel::pixToAng(double p, int s, int opticsGroup) const
{
	return s * angpix[opticsGroup] / p;
}

double ObservationModel::getPixelSize(int opticsGroup) const
{
	return angpix[opticsGroup];
}

int ObservationModel::numberOfOpticsGroups() const
{
	return opticsMdt.numberOfObjects();
}

bool ObservationModel::opticsGroupsSorted() const
{
	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		int og;
		opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);	
		
		if (og != i+1)
		{
			return false;
		}
	}
		
	return true;
}

std::vector<int> ObservationModel::findUndefinedOptGroups(const MetaDataTable &partMdt) const
{
	std::set<int> definedGroups;
	
	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		int og;
		opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);
		
		definedGroups.insert(og);
	}
	
	std::vector<int> out;
	out.reserve(opticsMdt.numberOfObjects());
	
	for (int i = 0; i < partMdt.numberOfObjects(); i++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);
		
		if (definedGroups.find(og) == definedGroups.end())
		{
			out.push_back(og);
		}
	}
	
	return out;
}

void ObservationModel::sortOpticsGroups(MetaDataTable& partMdt)
{
	std::map<int,int> old2new;
	
	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		int og;
		opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);
		
		old2new[og] = i+1;
		
		opticsMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, i+1, i);
	}
	
	for (int i = 0; i < partMdt.numberOfObjects(); i++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);
		partMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, old2new[og], i);
	}
}

std::vector<int> ObservationModel::getOptGroupsPresent(const MetaDataTable& partMdt) const
{
	const int gc = opticsMdt.numberOfObjects();
	const int pc = partMdt.numberOfObjects();
	
	std::vector<bool> optGroupIsPresent(gc, false);
	
	for (int p = 0; p < pc; p++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, p);
		
		optGroupIsPresent[og-1] = true;
	}
	
	std::vector<int> out(0);
	out.reserve(gc);
	
	for (int g = 0; g < gc; g++)
	{
		if (optGroupIsPresent[g])
		{
			out.push_back(g+1);
		}
	}
	
	return out;
}

const Image<Complex>& ObservationModel::getPhaseCorrection(int optGroup, int s)
{
	#pragma omp critical
	{
		if (phaseCorr[optGroup].find(s) == phaseCorr[optGroup].end())
		{
			if (phaseCorr[optGroup].size() > 100)
			{
				std::cerr << "Warning: " << (phaseCorr[optGroup].size()+1)
						  << " phase shift images in cache for the same ObservationModel." << std::endl;
			}
					
			const int sh = s/2 + 1;
			phaseCorr[optGroup][s] = Image<Complex>(sh,s);
			Image<Complex>& img = phaseCorr[optGroup][s];
			
			const double as = angpix[optGroup] * s;
			
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;
				
				for (int i = 0; i < oddZernikeCoeffs[optGroup].size(); i++)
				{
					int m, n;
					Zernike::oddIndexToMN(i, m, n);
					
					const double xx = x/as;
					const double yy = y < sh? y/as : (y-s)/as;
					
					phase += oddZernikeCoeffs[optGroup][i] * Zernike::Z_cart(m,n,xx,yy);
				}
				
				img(y,x).real = cos(phase);
				img(y,x).imag = sin(phase);
			}
		}
	}
	
	return phaseCorr[optGroup][s];
}

const Image<RFLOAT>& ObservationModel::getGammaOffset(int optGroup, int s)
{
	#pragma omp critical
	{
		if (gammaOffset[optGroup].find(s) == gammaOffset[optGroup].end())
		{
			if (gammaOffset[optGroup].size() > 100)
			{
				std::cerr << "Warning: " << (gammaOffset[optGroup].size()+1)
						  << " gamma offset images in cache for the same ObservationModel." << std::endl;
			}
					
			const int sh = s/2 + 1;
			gammaOffset[optGroup][s] = Image<RFLOAT>(sh,s);
			Image<RFLOAT>& img = gammaOffset[optGroup][s];
			
			const double as = angpix[optGroup] * s;
			
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;
				
				for (int i = 0; i < evenZernikeCoeffs[optGroup].size(); i++)
				{
					int m, n;
					Zernike::evenIndexToMN(i, m, n);
					
					const double xx = x/as;
					const double yy = y < sh? y/as : (y-s)/as;
					
					phase += evenZernikeCoeffs[optGroup][i] * Zernike::Z_cart(m,n,xx,yy);
				}
				
				img(y,x) = phase;
			}
		}	
	}
	
	return gammaOffset[optGroup][s];
}

bool ObservationModel::containsAllNeededColumns(const MetaDataTable& partMdt)
{
	return (partMdt.containsLabel(EMDL_ORIENT_ORIGIN_X)
         && partMdt.containsLabel(EMDL_ORIENT_ORIGIN_Y)
         && partMdt.containsLabel(EMDL_ORIENT_ROT)
         && partMdt.containsLabel(EMDL_ORIENT_TILT)
         && partMdt.containsLabel(EMDL_ORIENT_PSI)
         && partMdt.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET));
}
