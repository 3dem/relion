#ifndef FIT_BLOBS_PROGRAM_H
#define FIT_BLOBS_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/mesh/mesh_builder.h>

class Tomogram;
class TomogramSet;
class ManifoldSet;
class Blob3D;

class FitBlobs3DProgram
{
	public:

		
		FitBlobs3DProgram(){}

		
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



	private:


		void processTomogram(
				std::string tomoName,
				std::string spheresFn,
				TomogramSet& initial_tomogram_set,
				ManifoldSet& manifold_set);
		
		std::vector<double> segmentBlob(
				gravis::d3Vector sphere_position,
				double mean_radius_full,
				double radius_range,
				double binning,
				const RawImage<float>& preweighted_stack,
				double pixel_size,
				const std::vector<gravis::d4Matrix>& projections,
				const std::string& debug_prefix);
		
		Mesh createMesh(
				const std::vector<double>& blob_coeffs,
				double pixel_size,
				double spacing,
				double max_tilt_deg);

};

#endif
