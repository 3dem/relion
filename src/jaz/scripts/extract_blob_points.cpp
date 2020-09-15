#include <src/jaz/membrane/blob_fit_2d.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/util/drawing.h>
#include <src/args.h>
#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	std::string micrographs_list_filename, micrographs_dir, blobs_dir, outPath;
	double min_radius, max_radius, min_distance, perimeter_margin;
	int num_threads;
	bool diag, estimate_radius;
	
	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		micrographs_list_filename = parser.getOption("--ml", "Micrograph lists filename");
		micrographs_dir = parser.getOption("--md", "Micrographs directory");
		blobs_dir = parser.getOption("--bd", "Initial blobs directory");
		min_distance = textToDouble(parser.getOption("--d", "Min. particle distance [bin-1 pixels]", "128"));
		perimeter_margin = textToDouble(parser.getOption("--pm", "Distance of perimeter points from the outline [bin-1 pixels]", "64"));
		estimate_radius = parser.checkOption("--est_rad", "Estimate the radius of the blob");
		min_radius = textToDouble(parser.getOption("--r0", "Min. radius [bin-1 pixels]", "300"));
		max_radius = textToDouble(parser.getOption("--r1", "Max. radius [bin-1 pixels]", "800"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		
		diag = parser.checkOption("--diag", "Write out diagnostic information");
		
		outPath = parser.getOption("--o", "Output filename pattern");

		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			parser.writeUsage(std::cout);
			exit(1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	
	outPath = ZIO::makeOutputDir(outPath);
	
	const std::string perimeter_dir = "perimeter/Frames/";
	const std::string full_dir = "full/Frames/";
	const std::string outer_half_dir = "outer_half/Frames/";
	
	ZIO::makeOutputDir(outPath + perimeter_dir);
	ZIO::makeOutputDir(outPath + full_dir);
	ZIO::makeOutputDir(outPath + outer_half_dir);
	        
	micrographs_dir = ZIO::ensureEndingSlash(micrographs_dir);
	blobs_dir = ZIO::ensureEndingSlash(blobs_dir);
	
	std::vector<std::string> all_micrograph_names;
	all_micrograph_names.reserve(1000);
	
	std::ifstream file(micrographs_list_filename);
	
	if (!file)
	{
		REPORT_ERROR("Unable to read " + micrographs_list_filename);
	}
	
	std::string line;
	
	while (std::getline(file, line))
	{
		std::stringstream sts;
		sts << line;
		
		std::string name;
		sts >> name;
		
		all_micrograph_names.push_back(name);
	}
	
	const int micrograph_count = all_micrograph_names.size();
	
	std::vector<std::string> failed_micrographs;
		
	#pragma omp parallel for num_threads(num_threads)
	for (int m = 0; m < micrograph_count; m++)
	{
		const int thread_id = omp_get_thread_num();
		const bool verbose = thread_id == 0;
		
		if (verbose)
		{
			Log::beginSection(
				"Micrograph "+ZIO::itoa(m+1)+"/"
				+ZIO::itoa((int)std::ceil(micrograph_count / (double)num_threads))+" on this thread (out of "
				+ZIO::itoa(micrograph_count)+" on all threads)");
		}
		
		try 
		{
			const std::string micrograph_name = all_micrograph_names[m];
			const std::string micrograph_path = micrographs_dir + micrograph_name + ".mrc";
			const std::string blob_path = blobs_dir + micrograph_name + ".blobs";				
			std::vector<std::vector<double>> blob_shapes = ZIO::readDoublesTable(blob_path);
			
			BufferedImage<float> micrograph;			
			micrograph.read(micrograph_path);
			
			
			const double pixel_size = ImageFileHelper::getSamplingRate(micrograph_path);
			
			const int w = micrograph.xdim;
			const int h = micrograph.ydim;
			
			
			BufferedImage<float> weight(micrograph);
			weight.fill(1.f);
						
			const int blob_count = blob_shapes.size();
			
			BufferedImage<float> plot(micrograph);
			
			std::vector<DelineatedBlob2D> blobs(blob_count);
			
			if (verbose)
			{
				Log::beginProgress("Processing Blobs", blob_count);
			}
				            
			for (int b = 0; b < blob_count; b++)
			{
				if (verbose)
				{
					Log::updateProgress(b);
				}
				
				if (blob_shapes[b].size()%2 == 0)
				{
					for (int i = 0; i < blob_shapes[b].size(); i++)
					{
						std::cout << i << ": " << blob_shapes[b][i] << "\n";
					}
					
					REPORT_ERROR("The blobs in "+blob_path+" do not have a radius");
				}
				
				DelineatedBlob2D blob(blob_shapes[b]);
				
				if (estimate_radius)
				{				
					std::pair<std::vector<double>,std::vector<double>> radAvgAndWgh 
							= blob.radialAverageAndWeight(
								micrograph, weight, max_radius);
					
					std::vector<double> radAvg = radAvgAndWgh.first;
					std::vector<double> radWgh = radAvgAndWgh.second;
					
					const int bin_count = radAvg.size();
									
					int best_t = 0;
					double best_separation = 0.0;
					
					int delta_t = 20;
					
					for (int t = min_radius; t <= max_radius - delta_t && t < bin_count; t++)
					{
						double avg0 = 0.0;
						double wgh0 = 0.0;
						
						int t0 = t - delta_t;
						if (t0 < 0) t0 = 0;
						
						for (int i = t0; i < t; i++)
						{
							avg0 += radWgh[i] * radAvg[i];
							wgh0 += radWgh[i];
						}
						
						double avg1 = 0.0;
						double wgh1 = 0.0;
						
						int t1 = t + delta_t + 1;
						if (t0 > bin_count) t1 = bin_count;
						
						for (int i = t; i < t1; i++)
						{
							avg1 += radWgh[i] * radAvg[i];
							wgh1 += radWgh[i];
						}
						
						avg0 /= wgh0;
						avg1 /= wgh1;
						
						const double separation = avg1 - avg0;
						
						if (separation > best_separation)
						{
							best_separation = separation;
							best_t = t;
						}
					}
					
					blob.radius = best_t;
				}
				
				blobs[b] = blob;
			}
			
			if (verbose)
			{
				Log::endProgress();
			}
			
			const int EDGE = 0;
			const int OUTER_HALF = 1;
			const int FULL_BLOB = 2;
			
			std::vector<std::vector<d2Vector>> picks(3);
			
			
			for (int b = 0; b < blob_count; b++)
			{
				DelineatedBlob2D blob = blobs[b];
				
				
				const double d = min_distance;
				const double margin = d;
				const double dd = 0.5 * sqrt(3.0) * d;
				const int rc = (h - 2*margin) / dd + 1;
				const std::vector<int> cc = {(int)((w - 2*margin) / d), (int)((w - 2*margin - d/2) / d)};
				
				for (int r = 0; r < rc; r++)
				for (int c = 0; c < cc[r%2]; c++)
				{
					gravis::d2Vector p(
						margin + (r%2)*d/2 + c*d,
						margin + r*dd);
					
					double rr = blob.getRelativeSignedDistance(p);
					
					if (rr < 0)
					{
						picks[FULL_BLOB].push_back(p);
						
						if (diag) Drawing::drawCross(p, 100.f, 3, plot);
						
						if (rr > -0.5)
						{
							picks[OUTER_HALF].push_back(p);
							
							if (diag) Drawing::drawCross(p, 100.f, 6, plot);
						}
					}
				}
				
				const double l = blob.perimeter();
				const int samples = std::floor(l / min_distance);
				
				DelineatedBlob2D smaller_blob(blob);
				smaller_blob.radius -= perimeter_margin;
				
				for (int i = 0; i < samples; i++)
				{
					const double phi = 2 * PI * i / (double) samples;
					const d2Vector p = smaller_blob.getOutlinePoint(phi);
					
					if (p.x > margin && p.x < w - margin   
					 && p.y > margin && p.y < h - margin)
					{
						picks[EDGE].push_back(p);
						
						if (diag) Drawing::drawCross(p, 100.f, 10, plot);
					}
				}
			}
			
			if (diag) plot.write(outPath + "diag_"+ZIO::itoa(m)+".mrc", pixel_size);
			
			std::vector<MetaDataTable> tables(3);
			
			for (int i = 0; i < 3; i++)
			{
				const int pc = picks[i].size();
				MetaDataTable& table = tables[i];
				
				table.setName("images");
				
				table.addLabel(EMDL_IMAGE_COORD_X);
				table.addLabel(EMDL_IMAGE_COORD_Y);
				table.addLabel(EMDL_PARTICLE_AUTOPICK_FOM);
				table.addLabel(EMDL_PARTICLE_CLASS);
				table.addLabel(EMDL_ORIENT_PSI);
				
				for (int p = 0; p < pc; p++)
				{
					table.addObject();
					
					table.setValue(EMDL_IMAGE_COORD_X, picks[i][p].x, p);
					table.setValue(EMDL_IMAGE_COORD_Y, picks[i][p].y, p);
					table.setValue(EMDL_PARTICLE_AUTOPICK_FOM, 1.0, p);
					table.setValue(EMDL_PARTICLE_CLASS, 0, p);
					table.setValue(EMDL_ORIENT_PSI, 0.0, p);
				}
			}
			
			tables[EDGE].write(outPath + perimeter_dir + micrograph_name + "_blob_pick.star");
			tables[FULL_BLOB].write(outPath + full_dir + micrograph_name + "_blob_pick.star");
			tables[OUTER_HALF].write(outPath + outer_half_dir + micrograph_name + "_blob_pick.star");
		}
		catch (...)
		{
			#pragma omp critical
			{
				failed_micrographs.push_back(all_micrograph_names[m]);
			}
		}
		
		if (verbose)
		{
			Log::endSection();
		}
	}
	
	
	if (failed_micrographs.size() > 0)
	{
		ZIO::writeToFile(failed_micrographs, outPath + "failed_micrographs.txt");
	}
	
	return 0;
}
