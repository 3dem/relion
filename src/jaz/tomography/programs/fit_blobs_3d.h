#ifndef FIT_BLOBS_PROGRAM_H
#define FIT_BLOBS_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/tomography/optimisation_set.h>

class Tomogram;
class TomogramSet;
class ManifoldSet;
class Manifold;
class Blob3D;

class FitBlobs3DProgram
{
	public:

		
		FitBlobs3DProgram(){}

		
			std::string outPath;
			OptimisationSet optimisationSet;

			bool diag;

			int SH_bands, num_threads, max_iters;

			double
				relative_radius_range,
				membrane_separation,
				fiducials_radius_A,
				lowpass_sigma_real_A,
				fit_binning;



		void readParameters(int argc, char *argv[]);

		void run();



	private:


		void processTomogram(
				int tomo_index,
				const std::string& tomogram_name,
				const std::map<int, const Manifold*>& input_manifolds_map,
				const TomogramSet& tomogram_set,
				ManifoldSet& output_manifold_set);
		
		std::vector<double> segmentBlob(
				gravis::d3Vector sphere_position,
				double mean_radius_full,
				double radius_range,
				double membrane_separation,
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
