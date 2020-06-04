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
	IOParser parser;

	std::string in_phase, in_amplitude, in_blended, out_path;
	bool do_flip_phase;
	double filter_freq, output_scale;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		in_phase = parser.getOption("--phase", "Phase map");
		in_amplitude = parser.getOption("--amp", "Amplitude map");
		//in_blended = parser.getOption("--blend", "Consensus map");
		do_flip_phase = parser.checkOption("--flip_phase", "Flip the sign of the phase map");
		filter_freq = textToDouble(parser.getOption("--res", "Resolution [A]", "3.0"));
		output_scale = textToDouble(parser.getOption("--scale", "Scale of the output pixel colours", "0.02"));
		out_path = parser.getOption("--o", "Output path");

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	if (out_path[out_path.length()-1] != '/')
	{
		out_path = out_path + "/";
	}

	std::string command = "mkdir -p " + out_path;
	int res = system(command.c_str());

	Image<RFLOAT> dummy;
	dummy.read(in_phase, false);
	const double pixel_size = dummy.samplingRateX();

	std::cout << "pixel size: " << pixel_size << std::endl;

	std::vector<std::string> names {"blended", "phase", "amplitude"};

	const int blended_id = 0;
	const int phase_id = 1;
	const int amplitude_id = 2;

	std::vector<BufferedImage<double>> mapsRS(3);
	std::vector<BufferedImage<dComplex>> mapsFS(3);

	//mapsRS[blended_id].read(in_blended);
	mapsRS[phase_id].read(in_phase);
	mapsRS[amplitude_id].read(in_amplitude);

	if (do_flip_phase)
	{
		mapsRS[phase_id] *= -1.0;
	}

	for (int m = 1; m < 3; m++)
	{
		FFT::FourierTransform(mapsRS[m], mapsFS[m]);
	}


	const int s = mapsRS[1].xdim;
	const int sh = mapsFS[1].xdim;

	std::vector<std::vector<double>> power1D(3);


	for (int m = 1; m < 3; m++)
	{
		BufferedImage<double> power3D = PowerSpectrum::halfComplex(mapsFS[m]);

		power1D[m] = RadialAvg::fftwHalf_3D_lin(power3D);

		std::ofstream ofs(out_path + "power_" + names[m] + ".dat");

		for (int i = 0; i < power1D[m].size(); i++)
		{
			ofs << i << ' ' << power1D[m][i] << '\n';
		}
	}

	const int sc = power1D[1].size();

	std::vector<double> power_ratio(sc);

	{
		std::ofstream power_ratio_file(out_path + "power_ratio.dat");
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

	std::vector<double> optimal_scale(sc);

	{
		std::vector<double> num(sc, 0.0), denom(sc, 0.0);

		for (int z = 0; z < s;  z++)
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double r = RadialAvg::get1DIndex(x,y,z,s,s,s);
			const int ri = (int)(r+0.5);

			if (ri < sc)
			{
				const Complex za = mapsFS[amplitude_id](x,y,z);
				const Complex zp = mapsFS[phase_id](x,y,z);

				num[ri]   += za.real * zp.real + za.imag * zp.imag;
				denom[ri] += zp.real * zp.real + zp.imag * zp.imag;
			}
		}

		std::ofstream optimal_scale_file(out_path + "optimal_scale.dat");

		for (int i = 0; i < sc; i++)
		{
			if (denom[i] > 0.0)
			{
				optimal_scale[i] = num[i] / denom[i];
			}
			else
			{
				optimal_scale[i] = 0.0;
			}

			optimal_scale_file << i << ' ' << optimal_scale[i] << '\n';
		}
	}

	BufferedImage<dComplex> rescaled_amplitude_FS(sh,s,s);
	BufferedImage<dComplex> rescaled_phase_FS(sh,s,s);

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		double scale = RadialAvg::interpolate_FftwHalf_3D_lin(
						x,y,z, s,s,s, optimal_scale);

		rescaled_amplitude_FS(x,y,z) = mapsFS[amplitude_id](x,y,z) / scale;
		rescaled_phase_FS(x,y,z) = scale * mapsFS[phase_id](x,y,z);
	}

	BufferedImage<double> rescaled_phase_RS, rescaled_amplitude_RS;
	FFT::inverseFourierTransform(rescaled_phase_FS, rescaled_phase_RS);
	FFT::inverseFourierTransform(rescaled_amplitude_FS, rescaled_amplitude_RS);

	rescaled_phase_RS.write(out_path + "rescaled_phase.mrc", pixel_size);
	rescaled_amplitude_RS.write(out_path + "rescaled_amplitude.mrc", pixel_size);


	BufferedImage<double> difference_map = rescaled_phase_RS - mapsRS[amplitude_id];
	difference_map.write(out_path + "difference_map.mrc", pixel_size);


	const double resolution_Angstrom = filter_freq;
	const double resolution_pixels = (s * pixel_size) / resolution_Angstrom;


	BufferedImage<dComplex> filtered_phase_FS = ImageFilter::lowpass3D(
				mapsFS[phase_id], resolution_pixels, 10);

	BufferedImage<dComplex> filtered_amplitude_FS = ImageFilter::lowpass3D(
				mapsFS[amplitude_id], resolution_pixels, 10);

	BufferedImage<dComplex> filtered_rescaled_phase_FS = ImageFilter::lowpass3D(
				rescaled_phase_FS, resolution_pixels, 10);

	BufferedImage<dComplex> filtered_rescaled_amplitude_FS = ImageFilter::lowpass3D(
				rescaled_amplitude_FS, resolution_pixels, 10);


	BufferedImage<double>
			filtered_phase_RS,
			filtered_amplitude_RS,
			filtered_rescaled_phase_RS,
			filtered_rescaled_amplitude_RS;

	FFT::inverseFourierTransform(filtered_phase_FS, filtered_phase_RS);
	FFT::inverseFourierTransform(filtered_amplitude_FS, filtered_amplitude_RS);
	//FFT::inverseFourierTransform(filtered_rescaled_phase_FS, filtered_rescaled_phase_RS);
	//FFT::inverseFourierTransform(filtered_rescaled_amplitude_FS, filtered_rescaled_amplitude_RS);


	/*DualContrastWriter::writeAxialSlices(
		mapsRS[phase_id], mapsRS[amplitude_id], out_path + "dual_slice_", output_scale);

	DualContrastWriter::writeAxialSlices(
		mapsRS[phase_id], rescaled_amplitude_RS, out_path + "rescaled_dual_slice_", output_scale);

	DualContrastWriter::writeAxialSlices(
		filtered_phase_RS, filtered_amplitude_RS, out_path + "filtered_dual_slice_", output_scale);

	DualContrastWriter::writeAxialSlices(
		filtered_phase_RS, filtered_rescaled_amplitude_RS, out_path + "filtered_rescaled_dual_slice_", output_scale);*/

	filtered_phase_RS.write(out_path + "filtered_phase.mrc", pixel_size);
	filtered_amplitude_RS.write(out_path + "filtered_amplitude.mrc", pixel_size);


	const double high_pass_sigma = 20;

	BufferedImage<dComplex> high_pass_phase_FS = ImageFilter::highpassGauss3D(
				filtered_phase_FS, high_pass_sigma);

	BufferedImage<dComplex> high_pass_amplitude_FS = ImageFilter::highpassGauss3D(
				filtered_amplitude_FS, high_pass_sigma);

	BufferedImage<dComplex> high_pass_rescaled_phase_FS = ImageFilter::highpassGauss3D(
				filtered_rescaled_phase_FS, high_pass_sigma);


	FFT::inverseFourierTransform(high_pass_phase_FS, filtered_phase_RS);
	FFT::inverseFourierTransform(high_pass_amplitude_FS, filtered_amplitude_RS);
	FFT::inverseFourierTransform(high_pass_rescaled_phase_FS, filtered_rescaled_phase_RS);

	filtered_phase_RS.write(out_path + "high_pass_phase.mrc", pixel_size);
	filtered_amplitude_RS.write(out_path + "high_pass_amplitude.mrc", pixel_size);
	filtered_rescaled_phase_RS.write(out_path + "high_pass_rescaled_phase.mrc", pixel_size);

	(filtered_amplitude_RS - filtered_rescaled_phase_RS).write(
				out_path + "filtered_discrepancy.mrc", pixel_size);

	/*DualContrastWriter::writeAxialSlices(
		filtered_phase_RS, filtered_amplitude_RS, out_path + "high_pass_dual_slice_", output_scale);*/



	BufferedImage<double> correlation(sh,s,s), rot_correlation(sh,s,s);

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		dComplex za = mapsFS[amplitude_id](x,y,z);
		dComplex zp = mapsFS[phase_id](x,y,z);

		correlation(x,y,z) = za.real * zp.real + za.imag * zp.imag;
		rot_correlation(x,y,z) = za.real * zp.imag - za.imag * zp.real;
	}

	correlation.write(out_path + "correlation.mrc", pixel_size);
	rot_correlation.write(out_path + "rot_correlation.mrc", pixel_size);


	return 0;
}
