#ifndef MOTIF_DETECTION_H
#define MOTIF_DETECTION_H

#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/tomography/dynamo/catalogue.h>

class MotifDetection
{
	public:
		
		MotifDetection();
		
		MotifDetection(
			gravis::d4Matrix alignment,
			int centerIndex,
			int tomoIndex,
			double score);
		
		static Catalogue toCatalogue(
			const std::vector<MotifDetection>& occurences,
			double binningLevel);
		
		gravis::d4Matrix alignment;
		int centerIndex, tomoIndex;
		double score;
};

#endif
