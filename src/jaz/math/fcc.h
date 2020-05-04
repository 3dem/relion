#ifndef DYN_FCC_H
#define DYN_FCC_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/data_set.h>


class FCC
{
	public:
		
		static BufferedImage<double> compute(
				DataSet* dataSet,
				const std::vector<int>& partIndices,
				const Tomogram& tomogram,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				bool flip_value,
				int num_threads);
		
		static BufferedImage<double> compute3(
				DataSet* dataSet,
				const std::vector<int>& partIndices,
				const Tomogram& tomogram,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				bool flip_value,
				int num_threads);
		
		static BufferedImage<double> divide(
				const BufferedImage<double>& fcc3);
};

#endif
