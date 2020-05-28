#include <src/macros.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/color_helper.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_writer.h>
#include <src/jaz/math/fft.h>
#include <omp.h>


using namespace gravis;

int main(int argc, char *argv[])
{
	std::string in_path = "SpaBackproject/";
	std::string out_path = "dual_contrast/";
	std::string dc_root = "DC_SNR_0.003";
	std::string bc_root = "BC_SNR_0.003";

	const double pixel_size = 0.885;

	std::vector<std::string> names {"blended", "phase", "amplitude"};

	std::vector<std::string> suffix {"_merged.mrc", "_phase.mrc", "_amplitude.mrc"};
	std::vector<std::string> prefix {bc_root, dc_root, dc_root};



	const int blended_id = 0;
	const int phase_id = 1;
	const int amplitude_id = 2;


	std::vector<BufferedImage<double>> mapsRS(3);
	std::vector<BufferedImage<dComplex>> mapsFS(3);

	for (int m = 0; m < 3; m++)
	{
		mapsRS[m].read(in_path + prefix[m] + suffix[m]);
		FFT::FourierTransform(mapsRS[m], mapsFS[m]);
	}

	const int s = mapsRS[0].xdim;
	const int sh = mapsFS[0].xdim;

	std::vector<std::vector<double>> power1D(3);

	for (int m = 0; m < 3; m++)
	{
		BufferedImage<double> power3D = PowerSpectrum::halfComplex(mapsFS[m]);

		power1D[m] = RadialAvg::fftwHalf_3D_lin(power3D);

		std::ofstream ofs(out_path + "power_" + names[m] + ".dat");

		for (int i = 0; i < power1D.size(); i++)
		{
			ofs << i << ' ' << power1D[m][i] << '\n';
		}
	}

	const int sc = power1D[0].size();



	std::vector<double> power_ratio(sc);

	{
		std::ofstream power_ratio_file(out_path + "scale.dat");
		std::ofstream scale_file(out_path + "scale.dat");

		for (int i = 0; i < sc; i++)
		{
			const double phase_power = power1D[phase_id][i];
			const double amplitude_power = power1D[amplitude_id][i];

			power_ratio[i] = amplitude_power > 0.0?
						phase_power / amplitude_power : 0.0;

			power_ratio_file << i << ' ' << power_ratio[i] << '\n';
			scale_file << i << ' ' << sqrt(power_ratio[i]) << '\n';
		}
	}

	BufferedImage<dComplex> rescaled_amplitude_FS(sh,s,s);
	BufferedImage<dComplex> rescaled_phase_FS(sh,s,s);

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		double scale = sqrt(
					RadialAvg::interpolate_FftwHalf_3D_lin(
						x,y,z, s,s,s, power_ratio));

		rescaled_amplitude_FS(x,y,z) = scale * mapsFS[amplitude_id](x,y,z);
		rescaled_phase_FS(x,y,z) = mapsFS[phase_id](x,y,z) / scale;
	}

	BufferedImage<double> rescaled_phase_RS, rescaled_amplitude_RS;
	FFT::inverseFourierTransform(rescaled_phase_FS, rescaled_phase_RS);
	FFT::inverseFourierTransform(rescaled_amplitude_FS, rescaled_amplitude_RS);

	rescaled_phase_RS.write(out_path + "rescaled_phase.mrc", pixel_size);
	rescaled_amplitude_RS.write(out_path + "rescaled_amplitude.mrc", pixel_size);


	const double resolution_Angstrom = 3.0;
	const double resolution_pixels = (s * pixel_size) / resolution_Angstrom;


	BufferedImage<dComplex> filtered_phase_FS = ImageFilter::lowpass3D(
				mapsFS[phase_id], resolution_pixels, 10);

	BufferedImage<dComplex> filtered_amplitude_FS = ImageFilter::lowpass3D(
				mapsFS[amplitude_id], resolution_pixels, 10);

	BufferedImage<dComplex> filtered_rescaled_amplitude_FS = ImageFilter::lowpass3D(
				rescaled_amplitude_FS, resolution_pixels, 10);


	BufferedImage<double> filtered_phase_RS, filtered_amplitude_RS, filtered_rescaled_amplitude_RS;
	FFT::inverseFourierTransform(filtered_phase_FS, filtered_phase_RS);
	FFT::inverseFourierTransform(filtered_amplitude_FS, filtered_amplitude_RS);
	FFT::inverseFourierTransform(filtered_rescaled_amplitude_FS, filtered_rescaled_amplitude_RS);


	DualContrastWriter::writeAxialSlices(
		mapsRS[phase_id], mapsRS[amplitude_id], out_path + "dual_slice_", 0.1);

	DualContrastWriter::writeAxialSlices(
		mapsRS[phase_id], rescaled_amplitude_RS, out_path + "rescaled_dual_slice_", 0.1);

	DualContrastWriter::writeAxialSlices(
		filtered_phase_RS, filtered_amplitude_RS, out_path + "filtered_dual_slice_", 0.1);

	DualContrastWriter::writeAxialSlices(
		filtered_phase_RS, filtered_rescaled_amplitude_RS, out_path + "filtered_rescaled_dual_slice_", 0.1);

	return 0;
}
