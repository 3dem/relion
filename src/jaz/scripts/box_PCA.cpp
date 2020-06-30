#include <src/args.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/util/zio.h>
#include <src/Eigen/Dense>
#include <src/Eigen/SVD>

using namespace Eigen;
using namespace gravis;


int main(int argc, char *argv[])
{
	std::string all_input_filenames, output_directory;
	int number_of_dimensions_to_output;
	bool subtract_mean;
	
			
	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		all_input_filenames = parser.getOption(
			"--i", "Comma-separated list of input images (no space after comma)");
		
		number_of_dimensions_to_output = textToInteger(parser.getOption(
			"--N", "Number of dimensions to write out (negative means all)", "-1"));
		
		subtract_mean = parser.checkOption("--s", "Subtract the mean value from each image");
		
		output_directory = parser.getOption(
			"--o", "Output directory");
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	output_directory = ZIO::makeOutputDir(output_directory);
	
	std::vector<std::string> input_filenames = ZIO::split(all_input_filenames, ",");
	
	const int ic = input_filenames.size();	
	std::vector<BufferedImage<float>> images(ic);
	
	for (int i = 0; i < ic; i++)
	{
		std::cout << i << ": " << input_filenames[i] << std::endl;
		images[i].read(input_filenames[i]);
	}
	
	const int w = images[0].xdim;
	const int h = images[0].ydim;
	const int d = images[0].zdim;
	
	for (int i = 1; i < ic; i++)
	{
		if (!images[i].hasEqualSize(images[0]))
		{
			REPORT_ERROR_STR(
				"Images '" << input_filenames[0] << "' and '" << input_filenames[i] << "' are "
				<< "of different size: " << images[0].getSizeString() 
			    << " vs. " << images[i].getSizeString()  );
		}
	}
	
	
	std::vector<double> mean(ic);
	
	for (int i = 0; i < ic; i++)
	{
		mean[i] = subtract_mean? Normalization::computeMean(images[i]) : 0.0;
	}
	
	MatrixXd A0(ic,ic);
	
	std::vector<double> values(ic);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		for (int i = 0; i < ic; i++)
		{
			values[i] = images[i](x,y,z);
		}
		
		for (int i = 0; i < ic; i++)
		for (int j = 0; j < ic; j++)
		{
			A0(i,j) += (values[i] - mean[i]) * (values[j] - mean[j]);
		}
	}
	
	JacobiSVD<MatrixXd> svd(A0, ComputeFullV);
	
	const int N = number_of_dimensions_to_output > 0? number_of_dimensions_to_output : ic;

	for (int i = 0; i < N; i++)
	{
		BufferedImage<float> out = ((float)svd.matrixV()(0,i)) * images[0];
		
		for (int j = 1; j < ic; j++)
		{
			out += ((float)svd.matrixV()(j,i)) * images[j];
		}
		
		out.write(output_directory+ZIO::itoa(i)+".mrc");
	}
	
	return 0;
}
