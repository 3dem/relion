#ifndef BACKPROJECT_2D_H
#define BACKPROJECT_2D_H

#include <string>
#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_voxel.h>


class Backproject2D
{
	public:


		Backproject2D() { }


			std::string particlesFn, outDir;
			double margin, SNR;
			int num_threads;
			bool reextract, do_dual_contrast;


		void read(int argc, char **argv);
		void run();



	private:


		void backrotate_particle(
				const RawImage<fComplex>& image,
				long int particle_id,
				const MetaDataTable& particles_table,
				ObservationModel& obsModel,
				RawImage<dComplex> average,
				RawImage<double> weight);

		void backrotate_particle_dual_contrast(
				const RawImage<fComplex>& image,
				long int particle_id,
				const MetaDataTable& particles_table,
				ObservationModel& obsModel,
				RawImage<DualContrastVoxel<double>> average);


		BufferedImage<double> reconstruct(
				const RawImage<dComplex>& data,
				const RawImage<double>& weight,
				double Wiener_offset);

		std::pair<BufferedImage<double>, BufferedImage<double>>  reconstruct_dual_contrast(
				const RawImage<DualContrastVoxel<double>>& data,
				double Wiener_offset);


};

#endif
