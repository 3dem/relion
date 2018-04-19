#include <src/jaz/new_ft.h>
#include <src/image.h>

#include <omp.h>

void TakanoriTest()
{	
	int ny = 1000, nx = 1500;
	
	MultidimArray<double> real_array(ny, nx);
	MultidimArray<dComplex> complex_array;
	NewFFT::FourierTransform(real_array, complex_array);
	
	std::cout << real_array.data[5] << ", " << complex_array.data[5] << "\n";
}

void baseline()
{	
	int ny = 1000, nx = 1500;
	
	MultidimArray<double> real_array(ny, nx);
	MultidimArray<dComplex> complex_array;
	
	//NewFFT::FourierTransform(real_array, complex_array);
	FourierTransformer ft;
	ft.FourierTransform(real_array, complex_array);
	
	std::cout << real_array.data[5] << ", " << complex_array.data[5] << "\n";
}


int main(int argc, char *argv[])
{
	TakanoriTest();
	//baseline();
	return 0;
	
	Image<float> real(7420, 7676);
	
	#pragma omp parallel for num_threads(6)
	for (int i = 0; i < 100; i++)
	{
		Image<tComplex<float>> complex;
		
		std::cout << i << "\n";
		NewFFT::FourierTransform(real(), complex());
	}
	
	return 0;
}
