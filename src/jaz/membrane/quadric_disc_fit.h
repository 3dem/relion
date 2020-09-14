#ifndef MEMBRANE_DISC_FIT_H
#define MEMBRANE_DISC_FIT_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/optimization/optimization.h>

class QuadricDiscFit : public Optimization
{
	public:
		
		QuadricDiscFit(
			BufferedImage<float>& tiltSeries, 
			const std::vector<gravis::d4Matrix> proj, 
			gravis::d3Vector surface, 
			gravis::d3Vector inside,
			double bin,
			double diam,
			double thickness,
			int num_threads);

		BufferedImage<float> localTomo;
		gravis::d4Matrix localToTomo;
		
		
		
		double f(const std::vector<double> &x, void *tempStorage) const;
		
		std::vector<double> getInitial(double dist = 200.0);
		gravis::d4Matrix getMatrix(const std::vector<double> &x);
		
		BufferedImage<float> visualize(const std::vector<double> &x);
		
		BufferedImage<float> visualize2D(const std::vector<double> &x,
								 const BufferedImage<float>& stack,
								 const std::vector<gravis::d4Matrix>& proj);
		
		BufferedImage<float> writePhaselines(
				const std::vector<double> &x, int w, int h, 
				const std::vector<gravis::d4Matrix>& proj);
};

#endif
