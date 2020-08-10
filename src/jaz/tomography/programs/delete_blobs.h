#ifndef DELETE_BLOBS_PROGRAM_H
#define DELETE_BLOBS_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>

class Tomogram;
class Blob;

class DeleteBlobsProgram
{
	public:
		
		DeleteBlobsProgram(){}
		
			std::string outPath, tomoSetFn, tomoName, spheresFn, fiducialsDir;
			
			bool diag;
			int SH_bands, num_threads, max_iters;
			double sphere_thickness, spheres_binning, prior_sigma_A, fiducials_radius_A;

			std::vector<gravis::d4Vector> spheres;
			
		void readParameters(int argc, char *argv[]);
		void run();





	private:

		void processBlob(
				int blob_id,
				Tomogram& tomogram0,
				const std::vector<gravis::d3Vector>& fiducials);

		std::vector<gravis::d4Vector> readSpheresCMM(
				const std::string& filename,
				double binning);

		BufferedImage<float> drawTestFrame(
				Blob& blob,
				const Tomogram& tomogram,
				int test_frame,
				BufferedImage<float>& dummyWeight);
};

#endif
