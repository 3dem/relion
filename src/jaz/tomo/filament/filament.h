#ifndef FILAMENT_H
#define FILAMENT_H

#include <src/jaz/math/spline.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image/buffered_image.h>
#include <map>
#include <vector>

#include "filament_mapping.h"

class Filament
{
	public:
		
				static void computeArcLength(
						const Spline<gravis::d3Vector, double>& spline,
						int samples,
						std::map<double,double>& arc_to_index_out,
						std::map<double,double>& index_to_arc_out,
						double& length_out);
				
				static double translateIndex(
						double a, 
						const std::map<double,double>& table);
				
				static void computeSplineCoords(
						int w, int h,
						double binning,
						const Spline<gravis::d3Vector,double>& spline3D,
						const std::map<double,double>& index_to_arc,
						const std::map<double,double>& arc_to_index,
						double arcLen,
						const gravis::d4Matrix& proj,
						double maxDist,
						RawImage<float>& index_out,
						RawImage<float>& signed_dist_out,
						RawImage<float>& mask_out);
				
				static Spline<gravis::d2Vector,double> projectSpline(
						const Spline<gravis::d3Vector,double>& spline3D,
						const gravis::d4Matrix& proj);
				
				static gravis::d2Vector getSplineVelocity(
						const Spline<gravis::d2Vector,double>& spline,
						double t);
		
		
		
		
		Filament();		
		Filament(const Spline<gravis::d3Vector,double>& spline, double maxRadius);
		
		
			Spline<gravis::d3Vector,double> spline;
			double maxRadius;
			
			double arcLen;
			std::map<double,double> arc_to_index, index_to_arc;
			
			
		FilamentMapping rasteriseCoordinates(
				int w, int h, 
				const std::vector<gravis::d4Matrix>& proj,
				double binning,
				int num_threads) const;
			
			
};

#endif
