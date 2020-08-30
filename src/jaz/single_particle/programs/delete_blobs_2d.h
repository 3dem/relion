#ifndef DELETE_BLOBS_2D_PROGRAM_H
#define DELETE_BLOBS_2D_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>

class Blob2D;

class DeleteBlobs2DProgram
{
	public:

		
		DeleteBlobs2DProgram(){}
		
			std::string outPath, micrograph_filename, blobs_filename;
			
			bool diag;
			
			int max_frequencies, num_threads, max_iters;
			
			double 
				prior_sigma_A,
				highpass_sigma_real_A,
				max_binning, min_binning;

			std::vector<gravis::d4Vector> spheres;
			
		void readParameters(int argc, char *argv[]);

		void run();



	private:


		void processMicrograph(
				MetaDataTable& detected_blobs,
				const std::string& micrograph_name,
				BufferedImage<float>& visualisation,
				int micrograph_index,
				int micrograph_count);

		std::vector<double> fitBlob(
		        int blob_id,
				const std::vector<double>& initial_parameters,
		        double radius_full,
				double pixel_size_full,
				double binning_factor,
				const RawImage<float>& image_full,
				const std::string& image_name);

		BufferedImage<float> drawFit(
				Blob2D& blob,
				const BufferedImage<float>& image,
				const BufferedImage<float>& realWeight);

		BufferedImage<float> drawTestStack(
				Blob2D& blob,
				const BufferedImage<float>& image,
				const BufferedImage<float>& realWeight);

		std::vector<double> toBin1(
				const std::vector<double>& parameters,
				double binning_factor);

		std::vector<double> fromBin1(
				const std::vector<double>& parameters,
				double binning_factor);

		BufferedImage<float> evaluateRotationalSymmetry(
				const RawImage<float>& image,
				double radius, double max_radius, double sigma);
};

#endif
