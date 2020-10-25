#include "motif_detection.h"
#include <src/jaz/math/Euler_angles_dynamo.h>


using namespace gravis;

MotifDetection::MotifDetection()
{
}

MotifDetection::MotifDetection(d4Matrix alignment, int centerIndex, int tomoIndex, double score)
:	alignment(alignment),
	centerIndex(centerIndex),
	tomoIndex(tomoIndex),
	score(score)
{	
}

Catalogue MotifDetection::toCatalogue(
		const std::vector<MotifDetection>& occurences,
		double binningLevel)
{
	Catalogue out;
	
	const int pc = occurences.size();
	
	out.particles.resize(pc);
	
	for (int p = 0; p < pc; p++)
	{
		DynamoParticle& part = out.particles[p];
		const MotifDetection& d = occurences[p];
		
		part.x = binningLevel * d.alignment(0,3);
		part.y = binningLevel * d.alignment(1,3);
		part.z = binningLevel * d.alignment(2,3);
		
		part.dx = 0.0;
		part.dy = 0.0;
		part.dz = 0.0;
		
		part.tag = p;
		part.tomo = d.tomoIndex;
		
		d3Vector angs = EulerDynamo::matrixToAngles(d.alignment);
		
		part.tdrot = -RAD2DEG(angs[2]);
		part.tilt  = -RAD2DEG(angs[1]);
		part.narot = -RAD2DEG(angs[0]);
	}

	return out;	
}
