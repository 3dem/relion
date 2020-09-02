
#include <iostream>

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image/interpolation.h>


template <class T>
class Test
{
	public:
		
		Test();
		
		template <class T2>
		Test(const Test<T2>& z);
		
		T a, b;
};

template <class T>
Test<T>::Test()
{}

template <class T> template <class T2>
Test<T>::Test(const Test<T2>& z)
:	a(z.a), b(z.b)
{}





template<typename T> inline
T linearXY_clip(const RawImage<T>& img, double x, double y, int z)
{
	int x0 = FLOOR(x);
	int y0 = FLOOR(y);
	
	int x1 = x0 + 1;
	int y1 = y0 + 1;
	
	const double xf = x - x0;
	const double yf = y - y0;
	
	if (x0 < 0) x0 = 0;
	if (x0 >= img.xdim) x0 = img.xdim - 1;
	if (x1 < 0) x1 = 0;
	if (x1 >= img.xdim) x1 = img.xdim - 1;
	if (y0 < 0) y0 = 0;
	if (y0 >= img.ydim) y0 = img.ydim - 1;
	if (y1 < 0) y1 = 0;
	if (y1 >= img.ydim) y1 = img.ydim - 1;
	
	const T vx0 = (1 - xf) * img(x0,y0,z) + xf * img(x1,y0,z);
	const T vx1 = (1 - xf) * img(x0,y1,z) + xf * img(x1,y1,z);
	
	std::cout << vx0 << ", " << vx1 << std::endl;
	
	return (1 - yf) * vx0 + yf * vx1;
}






int main()
{
	const double x = 1.7, y = 2.3;
	const int z = 0;
	
	BufferedImage<tComplex<float>> img(5,5);
	img.fill(tComplex<float>(2,3));
	
	/*int x0 = FLOOR(x);
	int y0 = FLOOR(y);
	
	int x1 = x0 + 1;
	int y1 = y0 + 1;
	
	const double xf = x - x0;
	const double yf = y - y0;
	
	if (x0 < 0) x0 = 0;
	if (x0 >= img.xdim) x0 = img.xdim - 1;
	if (x1 < 0) x1 = 0;
	if (x1 >= img.xdim) x1 = img.xdim - 1;
	if (y0 < 0) y0 = 0;
	if (y0 >= img.ydim) y0 = img.ydim - 1;
	if (y1 < 0) y1 = 0;
	if (y1 >= img.ydim) y1 = img.ydim - 1;
	
	const tComplex<float> vx0 = (1 - xf) * img(x0,y0,z) + xf * img(x1,y0,z);
	const tComplex<float> vx1 = (1 - xf) * img(x0,y1,z) + xf * img(x1,y1,z);
	
	std::cout << (1 - yf) * vx0 + yf * vx1;*/
	
	std::cout << linearXY_clip(img.getRef(), x, y, z) << std::endl;
	
	return 0;
	
	Test<float> ti;
	ti.a = 5;
	ti.b = 7;
	
	Test<double> td = ti;
	
	std::cout << td.a << ", " << td.b << std::endl;
	
	return 0;
}
		
