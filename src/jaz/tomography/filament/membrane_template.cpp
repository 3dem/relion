#include "membrane_template.h"
#include <src/jaz/image/interpolation.h>


MembraneTemplate::MembraneTemplate(
		double spacing,  
		double inside_level, 
		double inside_halo_level,
		double membrane_level, 
		double centre_level,
		double outside_halo_level,
		double outside_level )

: d(spacing/2.0),
  spacing(spacing),
  inside_level(inside_level),
  inside_halo_level(inside_halo_level),
  membrane_level(membrane_level),
  centre_level(centre_level),
  outside_halo_level(outside_halo_level),
  outside_level(outside_level)
{	
}

double MembraneTemplate::getValue(double x)
{
	if (x < -3.0) // inside
	{
		return inside_level;
	}
	else if (x < -2.0)
	{
		return Interpolation::cosine(x + 3.0, inside_level, inside_halo_level);
	}
	else if (x < -1.0)
	{
		return Interpolation::cosine(x + 2.0, inside_halo_level, membrane_level);
	}
	else if (x < 0.0)
	{
		return Interpolation::cosine(x + 1.0, membrane_level, centre_level);
	}
	else if (x < 1.0)
	{
		return Interpolation::cosine(x,       centre_level, membrane_level);
	}
	else if (x < 2.0)
	{
		return Interpolation::cosine(x - 1.0, membrane_level, outside_halo_level);
	}
	else if (x < 3.0)
	{
		return Interpolation::cosine(x - 2.0, outside_halo_level, outside_level);
	}
	else
	{
		return outside_level;
	}
}
