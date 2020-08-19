#ifndef BACKPROJECT_2D_H
#define BACKPROJECT_2D_H

#include <string>
#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/single_particle/obs_model.h>


class Backproject2D
{
	public:


		Backproject2D() { }


			std::string particlesFn, outDir;


		void read(int argc, char **argv);
		void run();



	private:


		void backrotate_particle(const RawImage<fComplex> image,
				long int particle_id,
				const MetaDataTable& particles_table,
				ObservationModel &obsModel,
				RawImage<dComplex>& average,
				RawImage<double>& weight);

		BufferedImage<double> reconstruct(
				RawImage<dComplex>& data,
				RawImage<double>& weight,
				double Wiener_offset);


};

#endif
