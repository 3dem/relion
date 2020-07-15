#ifndef DUAL_CONTRAST_SOLUTION_H
#define DUAL_CONTRAST_SOLUTION_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t2Vector.h>

template <typename T>
class DualContrastSolution
{
	public:

		class ConditionInfo
		{
			public:

			ConditionInfo();

			T minimum, maximum, mean, std_deviation;
		};


		DualContrastSolution();
		DualContrastSolution(int w, int h, int d);


		BufferedImage<T> phase, amplitude;
		std::vector<ConditionInfo> conditionPerShell;
};

template <typename T>
DualContrastSolution<T>::ConditionInfo::ConditionInfo()
:	minimum(std::numeric_limits<T>::max()),
	maximum(-std::numeric_limits<T>::max()),
	mean(0),
	std_deviation(0)
{
}

template <typename T>
DualContrastSolution<T>::DualContrastSolution()
{
}

template <typename T>
DualContrastSolution<T>::DualContrastSolution(int w, int h, int d)
:	conditionPerShell(w/2+1)
{
}


#endif
