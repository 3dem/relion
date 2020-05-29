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
		in_blended = parser.getOption("--blend", "Consensus map");
		do_flip_phase = parser.checkOption("--flip_phase", "Flip the sign of the phase map");
		filter_freq = textToDouble(parser.getOption("--res", "Resolution [A]", "3.0"));
		output_scale = textToDouble(parser.getOption("--scale", "Scale of the output pixel colours", "0.03"));
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

	mapsRS[blended_id].read(in_blended);
	mapsRS[phase_id].read(in_phase);
	mapsRS[amplitude_id].read(in_amplitude);

	if (do_flip_phase)
	{
		mapsRS[phase_id] *= -1.0;
	}

	for (int m = 0; m < 3; m++)
	{
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


	const double resolution_Angstrom = filter_freq;
	const double resolution_pixels = (s * pixel_size) / resolution_Angstrom;


	BufferedImage<dComplex> filtered_phase_FS = ImageFilter::lowpass3D(
				mapsFS[phase_id], resolution_pixels, 10);

	BufferedImage<dComplex> filtered_amplitude_FS = ImageFilter::lowpass3D(
				mapsFS[amplitude_id], resolution_pixels, 10);

	BufferedImage<dComplex> filtered_rescaled_amplitude_FS = ImageFilter::lowpass3D(
				rescaled_amplitude_FS, resolution_pixels, 10);

	BufferedImage<dComplex> filtered_blended_FS = ImageFilter::lowpass3D(
				mapsFS[blended_id], resolution_pixels, 10);


	BufferedImage<double>
			filtered_phase_RS,
			filtered_amplitude_RS,
			filtered_rescaled_amplitude_RS,
			filtered_blended_RS;

	FFT::inverseFourierTransform(filtered_phase_FS, filtered_phase_RS);
	FFT::inverseFourierTransform(filtered_amplitude_FS, filtered_amplitude_RS);
	FFT::inverseFourierTransform(filtered_rescaled_amplitude_FS, filtered_rescaled_amplitude_RS);
	FFT::inverseFourierTransform(filtered_blended_FS, filtered_blended_RS);


	DualContrastWriter::writeAxialSlices(
		mapsRS[phase_id], mapsRS[amplitude_id], out_path + "dual_slice_", output_scale);

	DualContrastWriter::writeAxialSlices(
		mapsRS[phase_id], rescaled_amplitude_RS, out_path + "rescaled_dual_slice_", output_scale);

	DualContrastWriter::writeAxialSlices(
		filtered_phase_RS, filtered_amplitude_RS, out_path + "filtered_dual_slice_", output_scale);

	DualContrastWriter::writeAxialSlices(
		filtered_phase_RS, filtered_rescaled_amplitude_RS, out_path + "filtered_rescaled_dual_slice_", output_scale);


	d2Vector b(0.0, 0.0);
	d2Matrix A(0.0, 0.0, 0.0, 0.0);

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < s; x++)
	{
		const double ap = filtered_phase_RS(x,y,z);
		const double aa = filtered_rescaled_amplitude_RS(x,y,z);
		const double ab = filtered_blended_RS(x,y,z);

		// minimise |x_p * ap + x_a * aa - ab|²

		A(0,0) += ap * ap;
		A(0,1) += ap * aa;
		A(1,0) += aa * ap;
		A(1,1) += aa * aa;

		b[0] += ap * ab;
		b[1] += aa * ab;
	}

	d2Matrix Ainv = A;
	Ainv.invert();

	d2Vector optimum = Ainv * b;

	std::cout << "scales: " << optimum[0] << ", " << optimum[1]
			  << " (" << (100.0 * optimum[1] /
				 (std::abs(optimum[0]) + std::abs(optimum[1])))
			  << "%)" << std::endl;

	BufferedImage<double> recombined(s,s,s);
	BufferedImage<double> scaled_phase(s,s,s);
	BufferedImage<double> scaled_amplitude(s,s,s);

	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		recombined(x,y,z) =
			  optimum[0] * filtered_phase_RS(x,y,z)
			+ optimum[1] * filtered_rescaled_amplitude_RS(x,y,z);

		scaled_phase(x,y,z) =
			  optimum[0] * filtered_phase_RS(x,y,z);

		scaled_amplitude(x,y,z) =
			  optimum[1] * filtered_rescaled_amplitude_RS(x,y,z);
	}

	recombined.write(out_path + "recombined.mrc", pixel_size);
	filtered_blended_RS.write(out_path + "filtered_blended.mrc", pixel_size);
	scaled_phase.write(out_path + "filtered_phase.mrc", pixel_size);
	scaled_amplitude.write(out_path + "filtered_amplitude.mrc", pixel_size);

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		dComplex za = mapsFS[amplitude_id](x,y,z);
		dComplex zp = mapsFS[phase_id](x,y,z);


	}

	std::cout << std::endl;

	return 0;
}
