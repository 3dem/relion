#ifndef DYN_SPLINE_H
#define DYN_SPLINE_H

#include <vector>
#include <src/jaz/gravis/t4Vector.h>

template <typename T, typename S>
class Spline
{
	public:
		
		Spline();		
		Spline(std::vector<T> points);
		
			std::vector<T> 
				coeffs_a, 
				coeffs_b, 
				coeffs_c, 
				coeffs_d;
			
		T getValue(double t) const;
		T getTangent(double t) const;
		double length() const;
};

template <typename T, typename S>
Spline<T,S> :: Spline()
{}

template <typename T, typename S>
Spline<T,S> :: Spline(std::vector<T> points)
{
	const int len = points.size() - 1;
	
	coeffs_a.resize(len);
	coeffs_b.resize(len);
	coeffs_c.resize(len);
	coeffs_d.resize(len);
	
	for (int seg = 0; seg < len; seg++)
	{
		T p0 = points[seg];
		T p1 = points[seg+1];
		
		T v0 = (seg == 0)? p1 - p0 : (p1 - points[seg-1]) / S(2);
		T v1 = (seg == len-1)? p1 - p0 : (points[seg+2] - p1) / S(2);
		
		T u = p0 - p1 + v0;
		
		coeffs_a[seg] = v1 + S(2) * u;
		coeffs_b[seg] = -S(3) * u - v1;
		coeffs_c[seg] = v0;
		coeffs_d[seg] = p0;
	}
}
	
template <typename T, typename S>
T Spline<T,S> :: getValue(double t) const
{
	const int len = coeffs_a.size();
	
	int seg = (int) t;
	
	if (seg < 0) seg = 0;
	else if (seg >= len) seg = len - 1;
	
	const double dt = t - seg;
	
	const T a = coeffs_a[seg];
	const T b = coeffs_b[seg];
	const T c = coeffs_c[seg];
	const T d = coeffs_d[seg];
	
	return a * dt*dt*dt + b * dt*dt + c * dt + d;
}

template <typename T, typename S>
T Spline<T,S> :: getTangent(double t) const
{
	const int len = coeffs_a.size();
	
	int seg = (int) t;
	
	if (seg < 0) seg = 0;
	else if (seg >= len) seg = len - 1;
	
	const double dt = t - seg;
	
	const T a = coeffs_a[seg];
	const T b = coeffs_b[seg];
	const T c = coeffs_c[seg];
	const T d = coeffs_d[seg];
	
	return S(3) * a * dt*dt + S(2) * b * dt + c;
}

template <typename T, typename S>
double Spline<T,S> :: length() const
{
	return (double) coeffs_a.size();
}

#endif
