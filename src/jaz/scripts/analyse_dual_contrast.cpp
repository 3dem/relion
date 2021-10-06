#include <src/macros.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/color_helper.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/tomography/reconstruction.h>
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

double computeShellAverage(
		const RawImage<double>& map_RS,
		double r0, double r1)
{
	const int s = map_RS.xdim;

	double sum = 0.0;
	double count = 0.0;

	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const d3Vector r(
			x - s/2,
			y - s/2,
			z - s/2);

		const double rl = r.length();

		if (rl > r0 && rl < r1)
		{
			sum += map_RS(x,y,z);
			count += 1.0;
		}
	}

	return sum / count;
}

d2Vector computeUniformScaleAndIntercept(
		const RawImage<double>& phase_map_RS,
		const RawImage<double>& amp_map_RS)
{
	const int s = phase_map_RS.xdim;

	d2Matrix A(0,0,0,0);
	d2Vector b(0,0);

	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const d3Vector r(
			x - s/2,
			y - s/2,
			z - s/2);

		if (r.length() >= s/2 - 10) continue;

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

	return d2Vector(optimal_scale, optimal_offset);
}

double computeUniformScale(
		const RawImage<dComplex>& phase_map_FS,
		const RawImage<dComplex>& amp_map_FS)
{
	const int sh = phase_map_FS.xdim;
	const int s  = phase_map_FS.ydim;

	double num = 0.0, denom = 0.0;

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const d3Vector r(
			x - s/2,
			y - s/2,
			z - s/2);

		if (r.length() >= s/4) continue;

		const dComplex zp = phase_map_FS(x,y,z);
		const dComplex za = amp_map_FS(x,y,z);

		num   += za.imag * zp.imag + za.real * zp.real;
		denom += zp.imag * zp.imag + zp.real * zp.real;
	}

	return num / denom;
}

BufferedImage<dComplex> normalisePhaseByShell(
		const RawImage<dComplex>& phase_map_FS,
		const RawImage<dComplex>& amp_map_FS)
{
	const int sh = phase_map_FS.xdim;
	const int s  = phase_map_FS.ydim;

	std::vector<double> num(sh,0.0), denom(sh,0.0);

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const dComplex zp = phase_map_FS(x,y,z);
		const dComplex za = amp_map_FS(x,y,z);

		const double rd = RadialAvg::get1DIndex(x,y,z,s,s,s);
		const int r = ((int)(rd + 0.5) >= sh)? sh - 1 : (int)(rd + 0.5);

		num[r]   += za.imag * zp.imag + za.real * zp.real;
		denom[r] += zp.imag * zp.imag + zp.real * zp.real;
	}

	std::vector<double> ratio(sh);

	for (int r = 0; r < sh; r++)
	{
		if (denom[r] > 0.0)
		{
			ratio[r] = num[r] /= denom[r];
		}
		else
		{
			ratio[r] = 0.0;
		}
	}

	BufferedImage<dComplex> normalised_phase(sh,s,s);

	for (int z = 0; z < s;  z++)
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const double rd = RadialAvg::get1DIndex(x,y,z,s,s,s);
		const int r = ((int)(rd + 0.5) >= sh)? sh - 1 : (int)(rd + 0.5);

		normalised_phase(x,y,z) = ratio[r] * phase_map_FS(x,y,z);
	}

	return normalised_phase;
}

std::vector<d2Vector> extractCoords(
		const std::vector<d3Vector>& atom_positions,
		const RawImage<double> phase_map,
		const RawImage<double> amp_map,
		double coord_offset_pixels,
		double pixel_size,
		d2Matrix transform)
{
	const int ac = atom_positions.size();

	std::vector<d2Vector> out(ac);

	for (int a = 0; a < ac; a++)
	{

		const d3Vector pos = atom_positions[a] / pixel_size + coord_offset_pixels * d3Vector(1,1,1);

		const double vp = Interpolation::linearXYZ_clip(phase_map, pos.x, pos.y, pos.z);
		const double va = Interpolation::linearXYZ_clip(amp_map, pos.x, pos.y, pos.z);

		out[a] = transform * d2Vector(vp,va);
	}

	return out;
}

void plotCloud(
		CPlot2D& plot2D,
		const std::vector<d2Vector>& coords,
		int marker_size, dRGB colour, bool filled)
{
	CDataSet data_set;

	data_set.SetDrawMarker(true);
	data_set.SetDrawLine(false);
	data_set.SetMarkerSize(marker_size);
	data_set.SetDatasetColor(colour.r, colour.g, colour.b);
	data_set.SetDrawMarkerFilled(filled);

	const int cc = coords.size();

	for (int c = 0; c < cc; c++)
	{
		const d2Vector d = coords[c];
		data_set.AddDataPoint(CDataPoint(d.x, d.y));
	}

	plot2D.AddDataSet(data_set);
}

void plotEllipse(
		CPlot2D& plot2D,
		const Ellipse& ellipse, int samples,
		double line_thickness,
		int marker_size, dRGB colour, bool filled)
{
	CDataSet perimeter;

	perimeter.SetDrawMarker(false);
	perimeter.SetDrawLine(true);
	perimeter.SetDatasetColor(colour.r, colour.g, colour.b);
	perimeter.SetLineWidth(line_thickness);

	for (int i = 0; i < samples+1; i++)
	{
		const double phi = 2.0 * PI * i / (double) samples;
		const double sp = sin(phi);
		const double cp = cos(phi);

		const d2Vector c =
			ellipse.mean + 2.0 * (
				  cp * ellipse.axis0
				+ sp * ellipse.axis1);

		perimeter.AddDataPoint(CDataPoint(c.x, c.y));
	}

	plot2D.AddDataSet(perimeter);


	CDataSet centre;

	centre.SetDrawMarker(true);
	centre.SetDrawLine(false);
	centre.SetMarkerSize(marker_size);
	centre.SetDatasetColor(colour.r, colour.g, colour.b);
	centre.SetDrawMarkerFilled(filled);
	centre.AddDataPoint(CDataPoint(ellipse.mean.x, ellipse.mean.y));

	plot2D.AddDataSet(centre);


	CDataSet centre_outline;

	centre_outline.SetDrawMarker(true);
	centre_outline.SetDrawLine(false);
	centre_outline.SetMarkerSize(marker_size + 1);
	centre_outline.SetDatasetColor(0, 0, 0);
	centre_outline.SetDrawMarkerFilled(false);
	centre_outline.AddDataPoint(CDataPoint(ellipse.mean.x, ellipse.mean.y));

	plot2D.AddDataSet(centre_outline);
}

std::map<std::string, std::vector<d2Vector>> mergeElementClouds(
		const std::map<std::string, std::vector<d2Vector>>& name_to_cloud,
		const std::map<std::string, std::set<std::string>>& element_to_names)
{
	std::map<std::string, std::vector<d2Vector>> element_to_cloud;

	for (std::map<std::string, std::set<std::string>>::const_iterator it0 = element_to_names.begin();
		 it0 != element_to_names.end(); it0++)
	{
		const std::string element = it0->first;

		for (std::set<std::string>::const_iterator it1 = it0->second.begin();
			 it1 != it0->second.end(); it1++)
		{
			const std::vector<d2Vector>& coords0 = name_to_cloud.find(*it1)->second;
			std::vector<d2Vector>& coords = element_to_cloud[element];

			coords.insert(coords.begin(), coords0.begin(), coords0.end());
		}
	}

	return element_to_cloud;
}

Ellipse fitEllipse(const std::vector<d2Vector>& cloud)
{
	const int ac = cloud.size();

	d2Vector mean(0,0);

	for (int a = 0; a < ac; a++)
	{
		mean += cloud[a];
	}

	mean /= ac;

	Tensor2x2<double> cov(0,0,0);

	for (int a = 0; a < ac; a++)
	{
		const d2Vector d = cloud[a] - mean;

		cov.xx += d.x * d.x;
		cov.xy += d.x * d.y;
		cov.yy += d.y * d.y;
	}

	cov /= ac;

	return cov.getEllipse(mean);
}

double fitSlope(const std::vector<d2Vector>& cloud)
{
	const int ac = cloud.size();

	double num   = 0.0;
	double denom = 0.0;

	for (int a = 0; a < ac; a++)
	{
		num   += cloud[a].x * cloud[a].y;
		denom += cloud[a].x * cloud[a].x;
	}

	return num / denom;
}

void plotAllAtoms(
	std::map<std::string, std::set<std::string>> element_to_names,
	std::map<std::string,std::vector<d2Vector>> coords_by_name,
	std::map<std::string, dRGB> element_colours,
	double ellipse_line_width,
	bool plot_window_set,
	d2Vector plot_start,
	d2Vector plot_end,
	const std::string& out_path,
	const std::string& tag)
{
	std::ofstream element_angles_file(out_path+"element_angles_" + tag + ".txt");
	element_angles_file << "element  slope  angle [degrees]\n";

	CPlot2D all_atoms_plot("");

	std::map<std::string, std::vector<d2Vector>> all_element_coords;

	for (std::map<std::string, std::set<std::string>>::iterator it =
		 element_to_names.begin(); it != element_to_names.end(); it++)
	{
		std::string element = it->first;
		std::set<std::string>& names = it->second;

		std::vector<d2Vector>& all_coords = all_element_coords[element];

		for (std::set<std::string>::iterator it2 = names.begin();
			 it2 != names.end(); it2++)
		{
			const std::string atom_name = *it2;
			const std::vector<d2Vector>& coords = coords_by_name[atom_name];
			all_coords.insert(all_coords.end(), coords.begin(), coords.end());
		}

		dRGB full_colour = element_colours[element];
		const double fade = 0.67;
		dRGB faded_colour = (1 - fade) * full_colour + fade * dRGB(1,1,1);

		plotCloud(all_atoms_plot, all_coords, 1, faded_colour, true);

		const double slope = fitSlope(all_coords);
		const double angle_deg = RAD2DEG(atan(slope));

		element_angles_file << element << ": " << slope << ' ' << angle_deg << '\n';
	}

	for (std::map<std::string, std::vector<d2Vector>>::iterator it =
		 all_element_coords.begin(); it != all_element_coords.end(); it++)
	{
		std::string element = it->first;
		std::vector<d2Vector>& all_coords = it->second;

		Ellipse ellipse = fitEllipse(all_coords);

		dRGB colour = element_colours[element];

		CDataSet line;

		const double m = fitSlope(all_coords);
		const d2Vector optimal_ratio_spot = ellipse.mean.length() * d2Vector(1,m) / sqrt(1 + m*m);


		line.SetDrawMarker(false);
		line.SetDrawLine(true);
		line.SetDatasetColor(colour.r, colour.g, colour.b);
		line.AddDataPoint(CDataPoint(0, 0));
		line.AddDataPoint(CDataPoint(optimal_ratio_spot.x, optimal_ratio_spot.y));
		line.SetLineWidth(1.5);

		all_atoms_plot.AddDataSet(line);
	}

	for (std::map<std::string, std::vector<d2Vector>>::iterator it =
		 all_element_coords.begin(); it != all_element_coords.end(); it++)
	{
		std::string element = it->first;
		std::vector<d2Vector>& all_coords = it->second;

		plotEllipse(
			all_atoms_plot, fitEllipse(all_coords), 100,
			ellipse_line_width, 12, element_colours[element], true);
	}

	all_atoms_plot.SetTitle(tag);
	all_atoms_plot.SetXAxisTitle("phase");
	all_atoms_plot.SetYAxisTitle("amplitude");

	if (plot_window_set)
	{
		all_atoms_plot.SetViewArea(plot_start.x, plot_start.y, plot_end.x, plot_end.y);
	}

	all_atoms_plot.OutputPostScriptPlot(out_path+"all_atoms_" + tag + ".eps");
}



int main(int argc, char *argv[])
{
	IOParser parser;

	std::string in_phase, in_amplitude, in_model, out_path, tag;
	double filter_freq, high_pass_frequency, sigma_scale, psModel;
	d2Vector plot_start, plot_end;
	bool plot_window_set;
	int boxModel, boxOut, number_of_threads;
	bool normalise_by_shell, normalise_uniformly, write_filtered_maps, write_ratios, subtract_solvent;
	
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		in_phase = parser.getOption("--phase", "Phase map");
		in_amplitude = parser.getOption("--amp", "Amplitude map");
		in_model = parser.getOption("--pdb", "Atomic model");		
		filter_freq = textToDouble(parser.getOption("--res", "Resolution [A]", "3.0"));
		high_pass_frequency = textToDouble(parser.getOption("--hp", "High-pass frequency [px]", "-1.0"));
		sigma_scale = textToDouble(parser.getOption("--sc", "Region width for scale adaptation [px]", "7.0"));

		boxModel = textToInteger(parser.getOption("--box_model", "Box size of the map corresponding to the PDB file"));
		psModel = textToDouble(parser.getOption("--angpix_model", "Pixel size of the map corresponding to the PDB file"));
		boxOut = textToInteger(parser.getOption("--box_out", "Box size of the map to be compared"));
		
		subtract_solvent = parser.checkOption("--zp", "Zero solvent phase");
		normalise_uniformly = parser.checkOption("--nu", "Normalise phase map uniformly");
		normalise_by_shell = parser.checkOption("--ns", "Normalise phase map per shell");

		write_filtered_maps = parser.checkOption("--write_maps", "Write out filtered maps and their difference");
		write_ratios = parser.checkOption("--write_ratios", "Write out ratio maps");
		number_of_threads = textToInteger(parser.getOption("--j", "Number of threads", "6"));

		std::string plot_window = parser.getOption("--window", "Area to plot (format: <x0>,<x1>,<y0>,<y1>)", "");
		
		plot_window_set = plot_window != "";

		if (plot_window_set)
		{
			for (int i = 0; i < plot_window.length(); i++)
			{
				if (plot_window[i] == ',') plot_window[i] = ' ';
			}

			std::stringstream sts(plot_window);

			sts >> plot_start.x;
			sts >> plot_start.y;
			sts >> plot_end.x;
			sts >> plot_end.y;

			std::cout << "plotting from " << plot_start << " to " << plot_end << std::endl;
		}
		
		out_path = parser.getOption("--o", "Output path");
		tag = parser.getOption("--tag", "Output tag", "");

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


	const double resolution_pixels = filter_freq / pixel_size;
	
	phase_map_FS = ImageFilter::lowpass3D(phase_map_FS, resolution_pixels, 5);
	amp_map_FS = ImageFilter::lowpass3D(amp_map_FS, resolution_pixels, 5);
	
	if (high_pass_frequency > 0)
	{
		phase_map_FS = ImageFilter::highpassGauss3D(phase_map_FS, high_pass_frequency);
		amp_map_FS = ImageFilter::highpassGauss3D(amp_map_FS, high_pass_frequency);
	}

	FFT::inverseFourierTransform(phase_map_FS, phase_map_RS);
	FFT::inverseFourierTransform(amp_map_FS, amp_map_RS);

	if (subtract_solvent)
	{
		const double avg = computeShellAverage(phase_map_RS, s/4, s/2);
		phase_map_RS -= avg;
	}

	Reconstruction::taper(phase_map_RS, 10, false, number_of_threads);
	Reconstruction::taper(amp_map_RS, 10, false, number_of_threads);

	FFT::FourierTransform(phase_map_RS, phase_map_FS);
	FFT::FourierTransform(amp_map_RS, amp_map_FS);

	const double optimal_ratio = computeUniformScale(phase_map_FS, amp_map_FS);

	/*const d2Vector scale_and_intercept = computeUniformScaleAndIntercept(phase_map_RS, amp_map_RS);
	const double optimal_ratio = scale_and_intercept[0];
	const double optimal_intercept = scale_and_intercept[1];*/

	std::cout << "amplitude contrast: " << (100.0 * optimal_ratio) << '%' << std::endl;
	
	if (normalise_uniformly)
	{
		phase_map_FS *= optimal_ratio;
	}
	else if (normalise_by_shell)
	{
		phase_map_FS = normalisePhaseByShell(phase_map_FS, amp_map_FS);
	}
	
	
	FFT::inverseFourierTransform(phase_map_FS, phase_map_RS);
	FFT::inverseFourierTransform(amp_map_FS, amp_map_RS);
	
	if (write_filtered_maps)
	{
		BufferedImage<double> filtered_difference(s,s,s);

		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			filtered_difference(x,y,z) = amp_map_RS(x,y,z) - optimal_ratio * phase_map_RS(x,y,z);
		}


		phase_map_RS.write(out_path + "phase_" + tag + ".mrc", pixel_size);
		amp_map_RS.write(out_path + "amplitude_" + tag + ".mrc", pixel_size);
		filtered_difference.write(out_path + "difference_" + tag + ".mrc", pixel_size);
	}


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
		
		histogram.write(out_path + "joint_histogram_" + tag + ".mrc");
		
		BufferedImage<double> normalised_histogram(bins, bins);

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
				normalised_histogram(x,y) = histogram(x,y) / max_value;
			}
		}
		
		normalised_histogram.write(out_path + "joint_histogram_normalised_" + tag + ".mrc");

		BufferedImage<double> cumulative_histogram(bins, bins);

		for (int x = 0; x < bins; x++)
		{
			double sum = 0.0;

			for (int y = 0; y < bins; y++)
			{
				sum += histogram(x,y);
			}

			double sum_y = 0.0;

			for (int y = 0; y < bins; y++)
			{
				cumulative_histogram(x,y) = sum_y / sum;

				sum_y += histogram(x,y);
			}
		}

		cumulative_histogram.write(out_path + "cumulative_joint_histogram_" + tag + ".mrc");
	}

	const double coord_offset_pixels =
			( (boxOut/2)*pixel_size - (boxModel/2)*psModel ) / pixel_size;
	
	
	Assembly assembly;
	assembly.readPDB(in_model);

	std::set<std::string> heavy_external;
	heavy_external.insert("P");
	heavy_external.insert("CL");
	heavy_external.insert("FE");
	heavy_external.insert("ZN");
	heavy_external.insert("NA");
	heavy_external.insert("MG");

	std::map<std::string,std::vector<d3Vector>> atoms_by_name =
			PdbHelper::groupAtomsByName(assembly, heavy_external);

	std::map<std::string,std::vector<d2Vector>> coords_by_name;

	const d2Matrix scale(1e8, 0, 0, 1e8);

	for (std::map<std::string,std::vector<d3Vector>>::iterator it =
		 atoms_by_name.begin(); it != atoms_by_name.end(); it++)
	{
		coords_by_name[it->first] = extractCoords(
			it->second, phase_map_RS, amp_map_RS,
			coord_offset_pixels, pixel_size, scale);

		std::cout << it->first << ": " << coords_by_name[it->first].size() << std::endl;
	}

	std::map<std::string, std::set<std::string>> element_to_names;

	for (std::map<std::string,std::vector<d2Vector>>::iterator it =
		 coords_by_name.begin(); it != coords_by_name.end(); it++)
	{
		std::string atom_name = it->first;
		std::string element = PdbHelper::getElement(atom_name);

		element_to_names[element].insert(atom_name);
	}


	std::map<std::string, dRGB> element_colours;
	element_colours["C"] = dRGB(0.1,0.1,0.1);
	element_colours["N"] = dRGB(48,80,248)/255.0;
	element_colours["O"] = dRGB(1.0,0.0,0.0);
	element_colours["P"] = dRGB(0.5,0.125,1);
	element_colours["S"] = dRGB(0.9,0.9,0.1);
	element_colours["CL"] = dRGB(0.1,0.9,0.1);
	element_colours["FE"] = dRGB(255,165,0)/255.0;
	element_colours["ZN"] = dRGB(165,42,42)/255.0;
	element_colours["NA"] = dRGB(171,92,242)/255.0;
	element_colours["MG"] = dRGB(138,255,0)/255.0;


	const double ellipse_line_width = 0.5;

	{
		std::ofstream element_angles_file(out_path+"main_chain_element_angles_" + tag + ".txt");
		element_angles_file << "element  angle [degrees]\n";

		std::vector<d2Vector> C_coords  = coords_by_name[" C  "];
		std::vector<d2Vector> CA_coords = coords_by_name[" CA "];
		std::vector<d2Vector> O_coords  = coords_by_name[" O  "];
		std::vector<d2Vector> N_coords  = coords_by_name[" N  "];


		CPlot2D main_chain_plot("");

		std::vector<d2Vector> all_C_coords = C_coords;
		all_C_coords.insert(all_C_coords.end(), CA_coords.begin(), CA_coords.end());


		element_angles_file << "C: " << RAD2DEG(atan(fitSlope(all_C_coords))) << '\n';
		element_angles_file << "N: " << RAD2DEG(atan(fitSlope(N_coords))) << '\n';
		element_angles_file << "O: " << RAD2DEG(atan(fitSlope(O_coords))) << '\n';


		plotCloud(main_chain_plot, all_C_coords, 1, element_colours["C"], true);
		plotCloud(main_chain_plot, O_coords, 1, element_colours["O"], true);
		plotCloud(main_chain_plot, N_coords, 1, element_colours["N"], true);


		plotEllipse(main_chain_plot, fitEllipse(all_C_coords), 100,
					ellipse_line_width, 12, element_colours["C"], true);

		plotEllipse(main_chain_plot, fitEllipse(C_coords), 100,
					0.5*ellipse_line_width, 7, element_colours["C"], true);

		plotEllipse(main_chain_plot, fitEllipse(CA_coords), 100,
					0.5*ellipse_line_width, 7, element_colours["C"], true);

		plotEllipse(main_chain_plot, fitEllipse(O_coords), 100,
					ellipse_line_width, 12, element_colours["O"], true);

		plotEllipse(main_chain_plot, fitEllipse(N_coords), 100,
					ellipse_line_width, 12, element_colours["N"], true);

		main_chain_plot.SetXAxisTitle("phase");
		main_chain_plot.SetYAxisTitle("amplitude");

		if (plot_window_set)
		{
			main_chain_plot.SetViewArea(plot_start.x, plot_start.y, plot_end.x, plot_end.y);
		}

		main_chain_plot.OutputPostScriptPlot(out_path+"main_chain_atoms_" + tag + ".eps");
	}


	plotAllAtoms(
		element_to_names, coords_by_name, element_colours,
		ellipse_line_width, plot_window_set, plot_start, plot_end,
		out_path, tag);


	//plot2D.SetTitle("amplitude over phase");
	//plot2D.SetDrawLegend(true);

}
