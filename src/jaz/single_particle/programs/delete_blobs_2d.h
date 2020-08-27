#ifndef DELETE_BLOBS_2D_PROGRAM_H
#define DELETE_BLOBS_2D_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>

class Blob2D;

class DeleteBlobs2DProgram
{
	public:

		
		DeleteBlobs2DProgram(){}
		
			std::string outPath, micrographs_list_filename, blob_list_filename, fiducialsDir;
			
			bool diag;
			int max_frequencies, num_threads, max_iters, min_MG, max_MG;
			double blob_radius_A, blob_thickness_A,
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
				const std::vector<double>& initial,
				double outer_radius_full,
				double pixel_size_full,
				double binning_factor,
				RawImage<float>& image_full,
				int micrograph_index);

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
