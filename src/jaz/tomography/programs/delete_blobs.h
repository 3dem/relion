#ifndef DELETE_BLOBS_PROGRAM_H
#define DELETE_BLOBS_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>

class Tomogram;
class TomogramSet;
class Blob;

class DeleteBlobsProgram
{
	public:
		
		DeleteBlobsProgram(){}
		
			std::string outPath, tomoSetFn, listFn, fiducialsDir;
			
			bool diag;
			int SH_bands, num_threads, max_iters;
			double sphere_thickness, spheres_binning,
				prior_sigma_A, fiducials_radius_A,
				highpass_sigma_real_A,
				max_binning, min_binning;

			std::vector<gravis::d4Vector> spheres;
			
		void readParameters(int argc, char *argv[]);
		void run();
		
		void processTomogram(
		        std::string tomoName, 
		        std::string spheresFn,
		        TomogramSet& tomogramSet);




	private:

		std::vector<double> fitBlob(
				int blob_id,
				const std::vector<double>& initial,
				double binning_factor,
				Tomogram& tomogram0,
				const std::vector<gravis::d3Vector>& fiducials);

		std::vector<gravis::d4Vector> readSpheresCMM(
				const std::string& filename,
				double binning);

		BufferedImage<float> drawFit(
				Blob& blob,
				const Tomogram& tomogram,
				BufferedImage<float>& realWeight);

		BufferedImage<float> drawTestStack(
				Blob& blob,
				const Tomogram& tomogram,
				BufferedImage<float>& realWeight);
};

#endif
