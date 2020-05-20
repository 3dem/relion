#include "subtomo.h"
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/data_set.h>
#include <src/jaz/tomography/relion_data_set.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/optics/damage.h>
#include <src/time.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <iostream>

using namespace gravis;


void SubtomoProgram::run()
{
	DataSet* dataSet = DataSet::load(catFn, motFn);
	std::vector<std::vector<int> > particles = dataSet->splitByTomogram();
	
	if (cropSize < 0) cropSize = boxSize;
	
	const bool computeMultiplicity = write_multiplicity || write_normalised || write_combined;
	
	bool do_ctf = true;
	
	const long int tc = particles.size();
	const long int s2D = boxSize;
	const long int sh2D = s2D / 2 + 1;
	
	const long int s3D = cropSize;
	const long int sh3D = s3D / 2 + 1;
	
	const long int s02D = (int)(binning * s2D + 0.5);
	
	TomogramSet tomogramSet(tomoSetFn);

	if (dataSet->type == DataSet::Relion)
	{
		RelionDataSet* rds0 = (RelionDataSet*) dataSet;
		RelionDataSet rds = *rds0;
		
		for (int t = 0; t < tc; t++)
		{
			const int pc = particles[t].size();
			
			if (pc == 0) continue;
			
			for (int p = 0; p < pc; p++)
			{
				const int part_id = particles[t][p];
				
				const int opticsGroup = rds.getOpticsGroup(part_id);
				const double pixelSize = rds.getBinnedPixelSize(opticsGroup);
				
				std::string outData = outTag + "/" + dataSet->getName(part_id) + "_data.mrc";
				std::string outWeight = outTag + "/" + dataSet->getName(part_id) + "_weights.mrc";
				
				rds.setImageFileNames(outData, outWeight, part_id);
				
				d3Vector off, coord;
				
				rds.getParticleOffset(part_id, off.x, off.y, off.z);
				rds.getParticleCoord(part_id, coord.x, coord.y, coord.z);
				
				coord -= off / pixelSize;
				
				rds.setParticleOffset(part_id, 0,0,0);
				rds.setParticleCoord(part_id, coord.x, coord.y, coord.z);
			}
		}

		for (int og = 0; og < rds.numberOfOpticsGroups(); og++)
		{
			const double ps_img = rds.optTable.getDouble(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, og);
			const double ps_out = binning * ps_img;

			rds.optTable.setValue(EMDL_MICROGRAPH_BINNING, binning, og);
			rds.optTable.setValue(EMDL_MICROGRAPH_PIXEL_SIZE, ps_out, og);
			rds.optTable.setValue(EMDL_IMAGE_PIXEL_SIZE, ps_out, og);
			rds.optTable.setValue(EMDL_IMAGE_SIZE, cropSize, og);
		}

		rds.write(outTag + "_particles.star");
	}
		
	
	int dummy = std::system(("mkdir -p " + outTag).c_str());
	
	for (int t = 0; t < tc; t++)
	{
		const int pc = particles[t].size();		
		if (pc == 0) continue;		
		
		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));		
		Log::print("Loading");
				
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		
		const int fc = tomogram.frameCount;
		
		dataSet->checkTrajectoryLengths(particles[t][0], pc, fc, "subtomo");
		
		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s2D, binning);
		
		const int inner_thread_num = 1;
		const int outer_thread_num = num_threads / inner_thread_num;
		
		// @TODO: define input and output pixel sizes!
		
		const double binnedPixelSize = tomogram.optics.pixelSize * binning;
		
		Log::beginProgress("Backprojecting subtomograms", (int)ceil(pc/(double)outer_thread_num));
		
		#pragma omp parallel for num_threads(outer_thread_num)	
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();
			
			if (th == 0)
			{
				Log::updateProgress(p);
			}
						
			const int part_id = particles[t][p];
			
			const d3Vector pos = dataSet->getPosition(part_id);
			const std::vector<d3Vector> traj = dataSet->getTrajectoryInPix(
						part_id, fc, tomogram.optics.pixelSize);
			
			std::vector<d4Matrix> projCut(fc), projPart(fc);			
			
			BufferedImage<fComplex> particleStack = BufferedImage<fComplex>(sh2D,s2D,fc);			
			BufferedImage<float> weightStack(sh2D,s2D,fc);
			
			TomoExtraction::extractAt3D_Fourier(
					tomogram.stack, s02D, binning, tomogram.proj, traj,
					particleStack, projCut, inner_thread_num, false, true);
			
			
			if (!do_ctf) weightStack.fill(1.f);
			
			
			for (int f = 0; f < fc; f++)
			{
				if (do_rotate)
				{
					projPart[f] = projCut[f] * dataSet->getMatrix4x4(part_id, s3D, s3D, s3D);
				}
				else
				{
					projPart[f] = projCut[f];
				}
				
				if (do_ctf)
				{
					BufferedImage<float> ctfImg = TomoCtfHelper::drawCTF(
						tomogram.centralCTFs[f], tomogram.proj[f], pos, 
						tomogram.centre, tomogram.handedness, binnedPixelSize, s2D);
					
					const float scale = flip_value? -1.f : 1.f;
							
					for (int y = 0; y < s2D;  y++)
					for (int x = 0; x < sh2D; x++)
					{
						const float c = scale * ctfImg(x,y) * doseWeights(x,y,f);
						
						particleStack(x,y,f) *= c;
						weightStack(x,y,f) = c * c;
					}
				}
			}
			
			
			const int boundary = (boxSize - cropSize) / 2;
					
			if (boundary > 0)
			{
				BufferedImage<float> particlesRS, weightsRS;
				
				particlesRS = StackHelper::inverseFourierTransformStack(particleStack);
				weightsRS = StackHelper::inverseFourierTransformStack(FFT::toComplex(weightStack));
				
				particlesRS = Padding::unpadCenter2D_full(particlesRS, boundary);
				weightsRS = Padding::unpadCenter2D_full(weightsRS, boundary);
				
				TomoExtraction::cropCircle(particlesRS, 5, num_threads);
				TomoExtraction::cropCircle(weightsRS, 5, num_threads);
				
				particleStack = StackHelper::FourierTransformStack(particlesRS);
				weightStack = FFT::toReal(StackHelper::FourierTransformStack(weightsRS));
			}
			
			BufferedImage<fComplex> dataImgFS(sh3D,s3D,s3D);
			dataImgFS.fill(fComplex(0.0, 0.0));
			
			BufferedImage<float> ctfImgFS(sh3D,s3D,s3D), psfImgFS(sh3D,s3D,s3D), 
					dataImgRS(s3D,s3D,s3D), dataImgDivRS(s3D,s3D,s3D),
					multiImageFS;
			
			ctfImgFS.fill(0.0);
			psfImgFS.fill(0.0);
			dataImgRS.fill(0.0);			
			dataImgDivRS.fill(0.0);
			
			
			if (computeMultiplicity)
			{
				multiImageFS = BufferedImage<float>(sh3D,s3D,s3D);
				multiImageFS.fill(0.0);
				
				FourierBackprojection::backproject_bwd(
						particleStack, weightStack, projPart, dataImgFS,
						psfImgFS, ctfImgFS, multiImageFS, inner_thread_num);
			}
			else
			{
				FourierBackprojection::backproject_bwd(
						particleStack, weightStack, projPart, dataImgFS,
						psfImgFS, ctfImgFS, inner_thread_num);
			}
			
			Reconstruction::griddingCorrect3D(
					dataImgFS, psfImgFS, dataImgRS,
					true, inner_thread_num);
			
			
			if (do_cone_weight)
			{ 
				FFT::FourierTransform(dataImgRS, dataImgFS);
						
				d3Matrix R = dataSet->getMatrix3x3(part_id);
				
				for (int z = 0; z < s3D;  z++)
				for (int y = 0; y < s3D;  y++)
				for (int x = 0; x < sh3D; x++)
				{
					const d3Vector p0(
						x,
						y < s3D/2? y : y - s3D,
						z < s3D/2? z : z - s3D);
					
					const d3Vector p = R * p0;
					
					const double rho = sqrt(p.x*p.x + p.y*p.y);
					const double t = rho / (std::abs(p.z) * cone_slope + cone_sig0);
					
					const double m = 1.0 - exp(-0.5*t*t);
					
					dataImgFS(x,y,z) *= m;					
					psfImgFS(x,y,z) *= m;
					
					if (computeMultiplicity)
					{
						multiImageFS(x,y,z) *= m;
					}
				}
				
				FFT::inverseFourierTransform(dataImgFS, dataImgRS);
			}
				
			std::string outData = outTag + "/" + dataSet->getName(part_id) + "_data.mrc";
			std::string outWeight = outTag + "/" + dataSet->getName(part_id) + "_weights.mrc";
			std::string outCTF = outTag + "/" + dataSet->getName(part_id) + "_CTF2.mrc";
			std::string outDiv = outTag + "/" + dataSet->getName(part_id) + "_div.mrc";
			std::string outMulti = outTag + "/" + dataSet->getName(part_id) + "_multi.mrc";
			std::string outNrm = outTag + "/" + dataSet->getName(part_id) + "_data_nrm.mrc";
			std::string outWeightNrm = outTag + "/" + dataSet->getName(part_id) + "_CTF2_nrm.mrc";
						
			Reconstruction::taper(dataImgRS, taper, do_center, inner_thread_num);
			
			
			
			dataImgRS.write(outData);
			
			if (write_combined)
			{
				BufferedImage<float> ctfAndMultiplicity(sh3D,s3D,2*s3D);
				ctfAndMultiplicity.getSlabRef(0,s3D).copyFrom(ctfImgFS);
				ctfAndMultiplicity.getSlabRef(s3D,s3D).copyFrom(multiImageFS);
				
				ctfAndMultiplicity.write(outWeight);
			}
			
			if (write_ctf)
			{
				Centering::fftwHalfToHumanFull(ctfImgFS).write(outCTF);
			}
			
			if (write_multiplicity)
			{
				Centering::fftwHalfToHumanFull(multiImageFS).write(outMulti);
			}
			
			if (write_normalised)
			{
				BufferedImage<float> ctfImgFSnrm = ctfImgFS;
				BufferedImage<fComplex> dataImgCorrFS;	
				
				FFT::FourierTransform(dataImgRS, dataImgCorrFS, FFT::Both);
				
				for (long int i = 0; i < ctfImgFSnrm.getSize(); i++)
				{
					const float n = multiImageFS[i];
					ctfImgFSnrm[i] = n > 0.f? ctfImgFS[i] / n : 0.f;
					dataImgCorrFS[i] = n > 0.f? dataImgCorrFS[i] / n : fComplex(0.f,0.f);
				}
				
				FFT::inverseFourierTransform(dataImgCorrFS, dataImgDivRS, FFT::Both);
				
				dataImgDivRS.write(outNrm);
				Centering::fftwHalfToHumanFull(ctfImgFSnrm).write(outWeightNrm);
			}
			
			if (write_divided)
			{
				if (SNR > 0.0)
				{
					Reconstruction::ctfCorrect3D_Wiener(
						dataImgRS, ctfImgFS, dataImgDivRS,
						1.0 / SNR, inner_thread_num);
				}
				else
				{
					Reconstruction::ctfCorrect3D_heuristic(
						dataImgRS, ctfImgFS, dataImgDivRS,
						0.001, inner_thread_num);
				}
				
				Reconstruction::taper(dataImgDivRS, taper, do_center, inner_thread_num);				
				dataImgDivRS.write(outDiv);
			}
		}
		
		Log::endProgress();
		
		Log::endSection(); // tomogram
	}
}

BufferedImage<float> SubtomoProgram::cropAndTaper(const BufferedImage<float>& imgFS, int boundary, int num_threads) const
{
	BufferedImage<fComplex> ctfImgFS_complex = FFT::toComplex(imgFS);
	
	BufferedImage<float> ctfImgRS;
	FFT::inverseFourierTransform(ctfImgFS_complex, ctfImgRS);
	
	ctfImgRS = Centering::fftwFullToHumanFull(ctfImgRS);
	ctfImgRS = Padding::unpadCenter3D_full(ctfImgRS, boundary);
	
	Reconstruction::GaussEnvelope(ctfImgRS, env_sigma, do_center, num_threads);				
	Reconstruction::taper(ctfImgRS, taper, do_center, num_threads);
	
	BufferedImage<float> ctfImgRS_cent = Centering::humanFullToFftwFull(ctfImgRS);
	
	FFT::FourierTransform(ctfImgRS_cent, ctfImgFS_complex);
	
	return FFT::toReal(ctfImgFS_complex);
}
