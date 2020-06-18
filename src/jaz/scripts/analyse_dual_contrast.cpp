#include <src/macros.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/color_helper.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_writer.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/math/tensor2x2.h>
#include <src/jaz/atomic/pdb_helper.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/assembly.h>
#include <omp.h>
#include <set>


using namespace gravis;

void compareAtomsByResidue(
		const std::map<std::string,std::map<std::string,std::vector<d3Vector>>>& atoms_by_name_by_residue,
		double coord_offset_pixels,
        double pixel_size,
		const RawImage<double>& phase_map_RS,
		const RawImage<double>& amp_map_RS,
		std::string out_path,
		std::string tag,
		bool images_are_premultiplied);

void compareResidues(
		const std::map<std::string,std::vector<d3Vector>>& atoms_by_residue, 
		double coord_offset_pixels,
        double pixel_size,
		const RawImage<double>& phase_map_RS,
		const RawImage<double>& amp_map_RS,
		std::string filename_out,
		bool images_are_premultiplied);

int main(int argc, char *argv[])
{
	IOParser parser;

	std::string in_phase, in_amplitude, in_model, out_path;
	double filter_freq, high_pass_frequency, sigma_scale;
	
	double psModel;
	int boxModel, boxOut;	
	bool write_filtered_maps, write_ratios;
	
	
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
		
		write_filtered_maps = parser.checkOption("--write_maps", "Write out filtered maps and their difference");
		write_ratios = parser.checkOption("--write_ratios", "Write out ratio maps");
		
		
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
	
	//phase_map_FS *= optimal_ratio;
	
	
	FFT::inverseFourierTransform(phase_map_FS, phase_map_RS);
	FFT::inverseFourierTransform(amp_map_FS, amp_map_RS);
	
	
	{
		const int bins = 100;
		BufferedImage<double> histogram(bins, bins);
		histogram.fill(0.0);
		
		const double phase_sd = sqrt(Normalization::computeVariance(phase_map_RS, 0.0));
		const double amp_sd = sqrt(Normalization::computeVariance(amp_map_RS, 0.0));
		
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const double vp = phase_map_RS(x,y,z) / phase_sd;
			const double va = amp_map_RS(x,y,z) / amp_sd;
			
			const int ip = bins * (vp + 2.0) / 4.0 + 0.5;
			const int ia = bins * (va + 2.0) / 4.0 + 0.5;
			
			if (ip >= 0 && ip < bins && ia >= 0 && ia < bins)
			{
				histogram(ip,ia) += 1.0;
			}
		}
		
		histogram.write(out_path + "joint_histogram.mrc");
		
		for (int x = 0; x < bins; x++)
		{
			double max_value = 0.0;
			
			for (int y = 0; y < bins; y++)
			{
				if (histogram(x,y) > max_value) 
				{
					max_value = histogram(x,y);
				}
			}
			
			for (int y = 0; y < bins; y++)
			{
				histogram(x,y) /= max_value;
			}
		}
		
		histogram.write(out_path + "joint_histogram_normalised.mrc");
	}
	
	// find optimal scale and intercept in real space
	{
		d2Matrix A(0,0,0,0);
		d2Vector b(0,0);
		
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const double vp = phase_map_RS(x,y,z);
			const double va = amp_map_RS(x,y,z);
			
			A(0,0) += vp * vp;
			A(0,1) += vp *  1;
			A(1,0) +=  1 * vp;
			A(1,1) +=  1 *  1;
			
			b[0] += vp * va;
			b[1] +=      va;
		}
		
		d2Matrix Ainv = A;
		Ainv.invert();
		
		const d2Vector solution = Ainv * b;
		const double optimal_scale = solution[0];
		const double optimal_offset = solution[1];
		
		std::cout << "optimal scale = " << optimal_scale << std::endl;
		std::cout << "optimal offset = " << optimal_offset << std::endl;
	}
	
	const double coord_offset_A = (boxOut/2)*pixel_size - (boxModel/2)*psModel;
	const double coord_offset_pixels = coord_offset_A / pixel_size;
	
	
	Assembly assembly;
	assembly.readPDB(in_model);
	
	compareAtomsByResidue(
		PdbHelper::groupAtomsByNameByResidue(assembly),
		coord_offset_pixels, 
		pixel_size,
		phase_map_RS, amp_map_RS,
		out_path, "_sharp",
		false);
		
	compareResidues(
		PdbHelper::groupAtomsByResidue(assembly),
		coord_offset_pixels, 
		pixel_size,
		phase_map_RS, amp_map_RS,
		out_path + "residue_averages.txt",
		false);
		
	compareResidues(
		PdbHelper::groupAtomsByResidue(assembly, "C"),
		coord_offset_pixels, 
		pixel_size,
		phase_map_RS, amp_map_RS,
		out_path + "residue_averages_C.txt",
		false);
	
	return 0;
	
	if (write_filtered_maps)
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
			
			if (write_ratios) phase_2_smooth.write(out_path + "phase_2_smooth.mrc", pixel_size);	
			
			const double offset = 1e-16;
			
			if (write_ratios)		
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
			
			if (write_ratios)
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
			}
			
			if (write_ratios)
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
			
			const std::string tag = "_sigma_" + ZIO::itoa(sigma_scale);
			
			compareAtomsByResidue(
				PdbHelper::groupAtomsByNameByResidue(assembly),
				coord_offset_pixels, 
			    pixel_size,
				phase_2_smooth, amp_phase_smooth,
				out_path, tag,
				true);
				
			compareResidues(
				PdbHelper::groupAtomsByResidue(assembly),
				coord_offset_pixels, 
				pixel_size,
				phase_2_smooth, amp_phase_smooth,
				out_path + "residue_averages"+tag+".txt",
				true);
				
			compareResidues(
				PdbHelper::groupAtomsByResidue(assembly, "C"),
				coord_offset_pixels, 
				pixel_size,
				phase_2_smooth, amp_phase_smooth,
				out_path + "residue_averages_C"+tag+".txt",
				true);
			
		}
	}
		
	return 0;
}

void compareAtomsByResidue(
		const std::map<std::string,std::map<std::string,std::vector<d3Vector>>>& atoms_by_name_by_residue,
		double coord_offset_pixels,
		double pixel_size,
		const RawImage<double>& phase_map_RS,
		const RawImage<double>& amp_map_RS,
		std::string out_path,
		std::string tag,
		bool images_are_premultiplied)
{	
	std::ofstream averages(out_path + "atom_averages_by_residue"+tag+".txt");
	
	std::map<std::string, std::vector<d2Vector>> all_atoms;
	
	for (std::map<std::string,std::map<std::string,std::vector<d3Vector>>>::const_iterator itt = 
		 atoms_by_name_by_residue.begin(); itt != atoms_by_name_by_residue.end(); itt++)
	{
		const std::string residue_name = itt->first;
		const std::map<std::string,std::vector<d3Vector>>& atoms = itt->second;
		
		averages << residue_name << ":\n\n";
		averages << "atom   phase   ampl.   ratio   diff.\n\n";
		
		std::ofstream scatterplot, scatterplot_legend;
		
		if (!images_are_premultiplied)
		{
			scatterplot.open(out_path + "amp_over_phase" + tag + "_" + residue_name + ".dat");
			scatterplot_legend.open(out_path + "amp_over_phase" + tag + "_" + residue_name + "_legend.txt");
		}
		
		for (std::map<std::string,std::vector<d3Vector>>::const_iterator it = 
			 atoms.begin(); it != atoms.end(); it++)
		{
			const std::string atom_name = it->first;
			const std::vector<d3Vector>& positions = it->second;
			
			const int pc = positions.size();
			
			double num = 0.0;
			double denom = 0.0;
			double phase_sum = 0.0;
			double amp_sum = 0.0;
			
			for (int p = 0; p < pc; p++)
			{
				const d3Vector pos = positions[p] / pixel_size + coord_offset_pixels * d3Vector(1,1,1);
				
				const double vp = Interpolation::linearXYZ_clip(phase_map_RS, pos.x, pos.y, pos.z);
				const double va = Interpolation::linearXYZ_clip(amp_map_RS, pos.x, pos.y, pos.z);
				
				scatterplot << vp << ' ' << va << '\n';
				
				if (!images_are_premultiplied)
				{
					num += vp * va;
					denom += vp * vp;
					phase_sum += vp;
					amp_sum += va;
					
					scatterplot_legend << atom_name << '\n';
					scatterplot << '\n';
					
					all_atoms[atom_name].push_back(d2Vector(vp, va));
				}
				else
				{
					num += va;
					denom += vp;
				}
			}
			
			averages.setf( std::ios::fixed, std::ios::floatfield);
			
			if (!images_are_premultiplied)
			{
				averages 
					<< std::setw(3) << atom_name << "    " 
					<< std::setprecision(2) << std::setw(4) << phase_sum / (pc*1e-8) << "    " 
					<< std::setprecision(2) << std::setw(4) << amp_sum / (pc*1e-8) << "    " 
					<< std::setprecision(2) << std::setw(4) << num / denom << "    "  
					<< std::setprecision(2) << std::setw(4) << (amp_sum - phase_sum) / (pc*1e-8) << '\n'; 
			}
			else
			{
				averages 
					<< std::setw(3) << atom_name << " " 
					<< std::setprecision(2) << std::setw(4) << num / denom << '\n'; 
			} 
		}
		
		averages << "\n\n";
	}
	
	if (!images_are_premultiplied)
	{
		std::ofstream scatterplot_all(out_path + "amp_over_phase" + tag + "_all.dat");
		std::ofstream scatterplot_all_legend(out_path + "amp_over_phase" + tag + "_all_legend.txt");
		
		for (std::map<std::string,std::vector<d2Vector>>::iterator it = all_atoms.begin(); 
		     it != all_atoms.end(); it++)
		{
			const std::string atom_name = it->first;
			const std::vector<d2Vector>& coordinates = it->second;
			const int ac = coordinates.size();	
			
			scatterplot_all_legend << atom_name << " (" << ac << ")\n";
			
			d2Vector average(0,0);
			
			const bool rotate_coords = false;
			
			for (int a = 0; a < ac; a++)
			{
				const d2Vector c0 = coordinates[a];
				const d2Vector c = rotate_coords? d2Vector(c0.x+c0.y, c0.y-c0.x) : c0;
				        
				scatterplot_all << c.x << ' ' << c.y << '\n';
				
				average += c;
			}
			
			scatterplot_all << '\n';
			
			average /= ac;
			
			scatterplot_all << average.x << ' ' << average.y << "\n";
			scatterplot_all << "0 0\n\n";
						
			Tensor2x2<double> cov(0,0,0);
			
			for (int a = 0; a < ac; a++)
			{
				const d2Vector c00 = coordinates[a];
				const d2Vector c0 = rotate_coords? d2Vector(c00.x+c00.y, c00.y-c00.x) : c00;
				const d2Vector c = c0 - average;
				
				cov.xx += c[0] * c[0];
				cov.xy += c[0] * c[1];
				cov.yy += c[1] * c[1];
			}
			
			cov /= ac;
			
			d2Vector eigenvalues;
			d2Vector eigenvector0, eigenvector1;
			
						
			cov.diagonalize(eigenvalues, eigenvector0, eigenvector1);
			
			const d2Vector a0 = sqrt(eigenvalues[0]) * eigenvector0;
			const d2Vector a1 = sqrt(eigenvalues[1]) * eigenvector1;
				
			
			const int ellipse_samples = 99;
			
			for (int i = 0; i < ellipse_samples+1; i++)
			{
				const double phi = 2.0 * PI * i / (double) ellipse_samples;
				
				const d2Vector d0 = average + cos(phi) * a0 + sin(phi) * a1;
				
				scatterplot_all << d0.x << ' ' << d0.y << '\n';
			}
			
			scatterplot_all << '\n';
			
			for (int i = 0; i < ellipse_samples+1; i++)
			{
				const double phi = 2.0 * PI * i / (double) ellipse_samples;
				
				const d2Vector d0 = average + 2.0 * (cos(phi) * a0 + sin(phi) * a1);
				
				scatterplot_all << d0.x << ' ' << d0.y << '\n';
			}
			
			scatterplot_all << '\n';
		}
	}
}

void compareResidues(
		const std::map<std::string,std::vector<d3Vector>>& atoms_by_residue, 
		double coord_offset_pixels,
		double pixel_size,
		const RawImage<double>& phase_map_RS,
		const RawImage<double>& amp_map_RS,
		std::string filename_out,
		bool images_are_premultiplied)
{	
	std::ofstream averages(filename_out);
	
	averages << "res.   phase   ampl.   ratio   diff.\n\n";
	
	for (std::map<std::string,std::vector<d3Vector>>::const_iterator it = 
		 atoms_by_residue.begin(); it != atoms_by_residue.end(); it++)
	{
		const std::string residue_name = it->first;
		const std::vector<d3Vector>& positions = it->second;
					
		const int pc = positions.size();
		
		double num = 0.0;
		double denom = 0.0;
		double phase_sum = 0.0;
		double amp_sum = 0.0;
		
		for (int p = 0; p < pc; p++)
		{
			const d3Vector pos = positions[p] / pixel_size + coord_offset_pixels * d3Vector(1,1,1);
			
			const double vp = Interpolation::linearXYZ_clip(phase_map_RS, pos.x, pos.y, pos.z);
			const double va = Interpolation::linearXYZ_clip(amp_map_RS, pos.x, pos.y, pos.z);
			
			if (!images_are_premultiplied)
			{
				num += vp * va;
				denom += vp * vp;
				phase_sum += vp;
				amp_sum += va;
			}
			else
			{
				num += va;
				denom += vp;
			}
		}
		
		averages.setf( std::ios::fixed, std::ios::floatfield);
		
		if (!images_are_premultiplied)
		{
			averages 
				<< std::setw(3) << residue_name << "    " 
				<< std::setprecision(2) << std::setw(4) << phase_sum / (pc*1e-8) << "    " 
				<< std::setprecision(2) << std::setw(4) << amp_sum / (pc*1e-8) << "    " 
				<< std::setprecision(2) << std::setw(4) << num / denom << "    "  
				<< std::setprecision(2) << std::setw(4) << (amp_sum - phase_sum) / (pc*1e-8) << '\n'; 
		}
		else
		{
			averages 
				<< std::setw(3) << residue_name << " " 
				<< std::setprecision(2) << std::setw(4) << num / denom << '\n'; 
		}
	}
	
	averages << "\n";
}
	
