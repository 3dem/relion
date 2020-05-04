#ifndef FILAMENT_MEMBRANE_TEMPLATE_H
#define FILAMENT_MEMBRANE_TEMPLATE_H

class MembraneTemplate
{
	public:
		
		MembraneTemplate(
				double spacing,  // inter-membrane spacing (25-35A, ~40A here) 
				double inside_level, 
				double inside_halo_level,
				double membrane_level, 
				double centre_level,
				double outside_halo_level,
				double outside_level);
		
		
				const double 
						d,
						spacing,
						inside_level,
						inside_halo_level,
						membrane_level,
						centre_level,
						outside_halo_level,
						outside_level;
		
		
		double getValue(double x);
};

#endif
