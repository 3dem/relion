#include <src/jaz/new_ft.h>
#include <src/image.h>

int main(int argc, char *argv[])
{
	Image<float> real(7420, 7676);
	
	for (int i = 0; i < 100; i++)
	{
		Image<tComplex<float>> complex;
		
		std::cout << i << "\n";
		NewFFT::FourierTransform(real(), complex());
	}
	
	return 0;
}
