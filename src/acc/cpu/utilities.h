#ifndef CPU_UTILITIES_H
#define CPU_UTILITIES_H

namespace CpuKernels
{

template< typename T1, typename T2 >
static inline
int floorfracf(T1 a, T2 b)
{
	return (int)(a/b);
}

template< typename T1, typename T2 >
static inline
int ceilfracf(T1 a, T2 b)
{
	return (int)(a/b + 1);
}

} // end of namespace CpuKernels

#endif //CPU_UTILITIES_H
