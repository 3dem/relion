#include "spheroid.h"

using namespace gravis;


Spheroid::Spheroid(const std::vector<double> &parameters, int index)
: Manifold(index, getTypeName())
{
	for (int i = 0; i < 3; i++)
	{
		position[i] = parameters[i];
	}
	
	SH_coefficients.resize(parameters.size() - 3);
	
	for (int i = 3; i < parameters.size(); i++)
	{
		SH_coefficients[i-3] = parameters[i];
	}
}

std::vector<RigidAlignment> Spheroid::sampleParticles(
        double spacing, 
        double depth, 
        double minTilt, 
        double maxTilt)
{
	
}

std::vector<double> Spheroid::getParameters() const
{
	std::vector<double> parameters(SH_coefficients.size() + 3);
	        
	for (int i = 0; i < 3; i++)
	{
		parameters[i] = position[i];
	}
	
	for (int i = 3; i < parameters.size(); i++)
	{
		parameters[i] = SH_coefficients[i-3];
	}	
	
	return parameters;
}

std::string Spheroid::getTypeName()
{
	return "spheroid";
}
