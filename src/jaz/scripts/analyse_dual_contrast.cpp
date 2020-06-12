#include <src/macros.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/color_helper.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_writer.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/atomic/pdb_helper.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/assembly.h>
#include <omp.h>
#include <set>


using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;

	std::string in_phase, in_amplitude, in_model, out_path;
	double filter_freq, high_pass_frequency, sigma_scale;
	
	double psModel, psOut;
	int boxModel, boxOut;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		in_phase = parser.getOption("--phase", "Phase map");
		in_amplitude = parser.getOption("--amp", "Amplitude map");
		in_model = parser.getOption("--pdb", "Atomic model");		
		filter_freq = textToDouble(parser.getOption("--res", "Resolution [A]", "3.0"));
		high_pass_frequency = textToDouble(parser.getOption("--hp", "High-pass frequency [px]", "10.0"));
		sigma_scale = textToDouble(parser.getOption("--sc", "Region width for scale adaptation [px]", "7.0"));

		boxModel = textToInteger(parser.getOption("--box_model", "Box size of the map corresponding to the PDB file"));
		psModel = textToDouble(parser.getOption("--angpix_model", "Pixel size of the map corresponding to the PDB file"));
		boxOut = textToInteger(parser.getOption("--box_out", "Box size of the map to be compared"));
		
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


	BufferedImage<double> phase_map_RS(in_phase), amp_map_RS(in_amplitude);
	BufferedImage<dComplex> phase_map_FS, amp_map_FS;

	FFT::FourierTransform(phase_map_RS, phase_map_FS);
	FFT::FourierTransform(amp_map_RS, amp_map_FS);

	const int s = amp_map_RS.xdim;
	const int sh = amp_map_FS.xdim;


	const double resolution_pixels = (s * pixel_size) / filter_freq;
	
	phase_map_FS = ImageFilter::lowpass3D(phase_map_FS, resolution_pixels, 10);
	amp_map_FS = ImageFilter::lowpass3D(amp_map_FS, resolution_pixels, 10);
	
	phase_map_FS = ImageFilter::highpassGauss3D(phase_map_FS, high_pass_frequency);
	amp_map_FS = ImageFilter::highpassGauss3D(amp_map_FS, high_pass_frequency);
		
	
	double num = 0.0, denom = 0.0;
	
	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const dComplex zp = phase_map_FS(x,y,z);
		const dComplex za = amp_map_FS(x,y,z);
		
		num   += za.imag * zp.imag + za.real * zp.real;
		denom += zp.imag * zp.imag + zp.real * zp.real;
	}
	
	const double optimal_ratio = num / denom;
	
	std::cout << "amplitude contrast: " << (100.0 * optimal_ratio) << '%' << std::endl;
	
	phase_map_FS *= optimal_ratio;
	
	
	FFT::inverseFourierTransform(phase_map_FS, phase_map_RS);
	FFT::inverseFourierTransform(amp_map_FS, amp_map_RS);
	
	
	const double coord_offset = (boxOut/2)*pixel_size - (boxModel/2)*psModel;
	
	const bool merge_by_element = false;
	
	std::map<std::string,std::vector<d3Vector>> atoms = PdbHelper::splitAtoms(
				in_model, merge_by_element);
	
	{	
		std::ofstream averages(out_path + "element_averages.txt");
		
		averages << "element \tphase \tamplitude \topt. ratio\n";
		
		for (std::map<std::string,std::vector<d3Vector>>::iterator it = atoms.begin();
			 it != atoms.end(); it++)
		{
			const std::string element = it->first;
			const std::vector<d3Vector>& positions = it->second;
			
			std::ofstream scatterplot(out_path + "amp_over_phase_" + element + ".dat");
					
			const int pc = positions.size();
			
			double num = 0.0;
			double denom = 0.0;
			double phase_sum = 0.0;
			double amp_sum = 0.0;
			
			for (int p = 0; p < pc; p++)
			{
				const d3Vector pos = (positions[p] + coord_offset * d3Vector(1,1,1)) / pixel_size;
				
				const double vp = Interpolation::linearXYZ_clip(phase_map_RS, pos.x, pos.y, pos.z);
				const double va = Interpolation::linearXYZ_clip(amp_map_RS, pos.x, pos.y, pos.z);
				
				scatterplot << vp << ' ' << va << '\n';
				
				num += vp * va;
				denom += vp * vp;
				phase_sum += vp;
				amp_sum += va;
			}
			
			averages << element << " \t" 
					 << phase_sum / pc << " \t" 
					 << amp_sum / pc << " \t" 
					 << num / denom << '\n'; 
		}
	}
	
	return 0;
	
	{
		BufferedImage<double> filtered_difference = amp_map_RS - phase_map_RS;
		
		phase_map_RS.write(out_path + "phase.mrc", pixel_size);
		amp_map_RS.write(out_path + "amplitude.mrc", pixel_size);
		filtered_difference.write(out_path + "difference.mrc", pixel_size);
		
		const double phase_variance = Normalization::computeVariance(phase_map_RS, 0.0);
		const double amp_variance = Normalization::computeVariance(amp_map_RS, 0.0);
		const double diff_variance = Normalization::computeVariance(amp_map_RS, 0.0);
		
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const double a = amp_map_RS(x,y,z);
			const double p = phase_map_RS(x,y,z);
			const double d = filtered_difference(x,y,z);
			
			const double ma = a > 0.0? 1.0 - exp( -0.5 * a*a / amp_variance ) : 0.0;
			const double mp = p > 0.0? 1.0 - exp( -0.5 * p*p / phase_variance ) : 0.0;
			const double md = d > 0.0? 1.0 - exp( -0.5 * d*d / diff_variance ) : 0.0;
			
			filtered_difference(x,y,z) *= ma * mp * md;
		}
		
		filtered_difference.write(out_path + "masked_difference.mrc", pixel_size);
	}
	
	
	{
		BufferedImage<double> amp_phase_RS(s,s,s), phase_2_RS(s,s,s);
		
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			amp_phase_RS(x,y,z) = amp_map_RS(x,y,z) * phase_map_RS(x,y,z);
			phase_2_RS(x,y,z) = phase_map_RS(x,y,z) * phase_map_RS(x,y,z);
		}
		
		const bool estimate_intercept = false;
		
		if (estimate_intercept)
		{
			BufferedImage<double> amp_smooth = ImageFilter::Gauss3D(amp_map_RS, sigma_scale);
			BufferedImage<double> phase_smooth = ImageFilter::Gauss3D(amp_map_RS, sigma_scale);
			BufferedImage<double> amp_phase_smooth = ImageFilter::Gauss3D(amp_phase_RS, sigma_scale);
			BufferedImage<double> phase_2_smooth = ImageFilter::Gauss3D(phase_2_RS, sigma_scale);
			
			BufferedImage<double> local_scale(s,s,s), local_intercept(s,s,s);
					
			for (int z = 0; z < s; z++)
			for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				const double a = amp_smooth(x,y,z);
				const double p = phase_smooth(x,y,z);
				const double ap = amp_phase_smooth(x,y,z);
				const double p2 = phase_2_smooth(x,y,z);
				
				d2Matrix A(p2,p,p,1);
				const d2Vector b(ap,a);
				
				A.invert();
				
				const d2Vector opt = A * b;
				
				local_scale(x,y,z) = opt[0];
				local_intercept(x,y,z) = opt[1];
			}
			
			local_scale.write(out_path + "local_scale.mrc", pixel_size);
			local_intercept.write(out_path + "local_intercept.mrc", pixel_size);
			
			BufferedImage<double> local_difference(s,s,s);
			
			for (int z = 0; z < s; z++)
			for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				const double a = amp_map_RS(x,y,z);
				const double p = phase_map_RS(x,y,z);
				const double m = local_scale(x,y,z);
				const double q = local_intercept(x,y,z);
				
				local_difference(x,y,z) = a - (m*p + q);
			}
			
			local_difference.write(out_path + "local_difference.mrc", pixel_size);
		}
		else
		{
			BufferedImage<double> amp_phase_smooth = ImageFilter::Gauss3D(amp_phase_RS, sigma_scale);
			BufferedImage<double> phase_2_smooth = ImageFilter::Gauss3D(phase_2_RS, sigma_scale);
			
			phase_2_smooth.write(out_path + "phase_2_smooth.mrc", pixel_size);	
			
			const double offset = 1e-16;
					
			{
				BufferedImage<double> local_difference(s,s,s);
				
				for (int z = 0; z < s; z++)
				for (int y = 0; y < s; y++)
				for (int x = 0; x < s; x++)
				{
					const double a = amp_map_RS(x,y,z);
					const double p = phase_map_RS(x,y,z);
					const double ap = amp_phase_smooth(x,y,z);
					const double p2 = phase_2_smooth(x,y,z) + offset;
					const double scale = p2 > 0.0? ap / p2 : 0.0;
					
					local_difference(x,y,z) = a - scale * p;
				}
				
				local_difference.write(out_path + "local_difference.mrc", pixel_size);			
			}
			
			{
				BufferedImage<double> local_scale(s,s,s);
				
				for (int z = 0; z < s; z++)
				for (int y = 0; y < s; y++)
				for (int x = 0; x < s; x++)
				{
					const double ap = amp_phase_smooth(x,y,z);
					const double p2 = phase_2_smooth(x,y,z) + offset;
					const double scale = p2 > 0.0? ap / p2 : 0.0;
					
					local_scale(x,y,z) = scale;
				}
				
				local_scale.write(out_path + "local_scale.mrc", pixel_size);
				
				for (std::map<std::string,std::vector<d3Vector>>::iterator it = atoms.begin();
					 it != atoms.end(); it++)
				{
					const std::string element = it->first;
					const std::vector<d3Vector>& positions = it->second;
					
					std::ofstream local_scale_file(out_path + "local_scale_" + element + ".dat");
							
					const int pc = positions.size();
					
					for (int p = 0; p < pc; p++)
					{
						const d3Vector pos = (positions[p] + coord_offset * d3Vector(1,1,1)) / pixel_size;
						
						const double vsp = Interpolation::linearXYZ_clip(local_scale, pos.x, pos.y, pos.z);
						
						local_scale_file << vsp << " 0\n";
					}
				}		
			}
			
			{
				BufferedImage<double> scaled_phase(s,s,s);
				
				for (int z = 0; z < s; z++)
				for (int y = 0; y < s; y++)
				for (int x = 0; x < s; x++)
				{
					const double p = phase_map_RS(x,y,z);
					const double ap = amp_phase_smooth(x,y,z);
					const double p2 = phase_2_smooth(x,y,z) + offset;
					const double scale = p2 > 0.0? ap / p2 : 0.0;
					
					scaled_phase(x,y,z) = scale * p;
				}
				
				scaled_phase.write(out_path + "scaled_phase.mrc", pixel_size);
			}
			
			const double coord_offset = (boxOut/2)*pixel_size - (boxModel/2)*psModel;
			
			std::map<std::string,std::vector<d3Vector>> atoms = PdbHelper::splitAtoms(
						in_model, merge_by_element);
			
			{	
				std::ofstream averages(out_path + "element_averages_regional.txt");
				
				averages << "element \topt. ratio\n";
				
				for (std::map<std::string,std::vector<d3Vector>>::iterator it = atoms.begin();
					 it != atoms.end(); it++)
				{
					const std::string element = it->first;
					const std::vector<d3Vector>& positions = it->second;
					
					const int pc = positions.size();
					
					double num = 0.0;
					double denom = 0.0;
					
					for (int p = 0; p < pc; p++)
					{
						const d3Vector pos = (positions[p] + coord_offset * d3Vector(1,1,1)) / pixel_size;
						
						const double vn = Interpolation::linearXYZ_clip(amp_phase_smooth, pos.x, pos.y, pos.z);
						const double vd = Interpolation::linearXYZ_clip(phase_2_smooth, pos.x, pos.y, pos.z);
						
						num += vn;
						denom += vd;
					}
					
					averages << element << " \t" 
							 << num / denom << '\n'; 
				}
			}
			
		}
	}
		
	return 0;
}
