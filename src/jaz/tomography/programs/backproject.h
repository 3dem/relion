#ifndef BACKPROJECT_PROGRAM_H
#define BACKPROJECT_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>


class BackprojectProgram
{
	public:
		
		BackprojectProgram(){}
		
			std::string outTag, tomoSetFn, catFn, motFn, symmName;
			
			bool do_whiten, diag, no_subpix_off, explicit_gridding, no_reconstruction;
			int boxSize, cropSize, num_threads, outer_threads, inner_threads, max_mem_GB;
			double SNR, taper, binning;
			
		void readParameters(int argc, char *argv[]);
		void run();


	protected:

		void reconstruct(
				BufferedImage<double>& dataImgRS,
				BufferedImage<double>& dataImgDivRS,
				BufferedImage<double>& ctfImgFS,
				BufferedImage<double>* psfImgFS,
				BufferedImage<dComplex>& dataImgFS);

		void writeOutput(
				const BufferedImage<double>& corrected,
				const BufferedImage<double>& data,
				const BufferedImage<double>& weight,
				const std::string& tag,
				double pixelSize);
};

#endif
