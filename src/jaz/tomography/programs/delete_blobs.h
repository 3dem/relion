#ifndef DELETE_BLOBS_PROGRAM_H
#define DELETE_BLOBS_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>


class DeleteBlobsProgram
{
	public:
		
		DeleteBlobsProgram(){}
		
			std::string outPath, tomoSetFn, tomoName, spheresFn;
			
			bool diag;
			int SH_bands, num_threads, max_iters;
			double inner_margin, outer_margin, spheres_binning;

			std::vector<gravis::d4Vector> spheres;
			
		void readParameters(int argc, char *argv[]);
		void run();





	private:

		std::vector<gravis::d4Vector> readSpheresCMM(
				const std::string& filename,
				double binning);
};

#endif
