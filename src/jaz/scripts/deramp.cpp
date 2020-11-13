#include <src/jaz/math/Euler_angles_dynamo.h>
#include <src/macros.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/math/fft.h>
#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{
    std::vector<std::string> indices = {"R2_3115", "R2_3160", "R2_3180", "R2_3197"};

    for (int ind = 0; ind < indices.size(); ind++)
    {
        std::string index = indices[ind];

        std::string fn = index+"_square_DW_4.1_CTFcorrected_FlippedXZY_WeightFiltered_XYZ_bin8.rec:mrc";
        std::string fn_out = index+"_square_DW_4.1_CTFcorrected_FlippedXZY_WeightFiltered_XYZ_bin8_filtered.rec:mrc";
        std::string fn_test = index+"_square_DW_4.1_CTFcorrected_FlippedXZY_WeightFiltered_XYZ_bin8_unfiltered.rec:mrc";

        std::cout << "   " << fn << "\n->   " << fn_out << std::endl;
        const double r_crit = 50;

        BufferedImage<float> img;
        img.read(fn);

        img.write(fn_test);

        BufferedImage<fComplex> imgFS;

        FFT::FourierTransform(img, imgFS, FFT::Both);

        const int w = img.xdim;
        const int h = img.ydim;
        const int d = img.zdim;

        const int wh = w/2 + 1;

        for (int z = 0; z < d; z++)
        for (int y = 0; y < h; y++)
        for (int x = 0; x < wh; x++)
        {
            const double xx = x;
            const double yy = y < h/2? y : y - h;
            const double zz = z < d/2? z : z - d;

            const double r = sqrt(xx*xx + yy*yy + zz*zz);

            const double amp = r_crit / (r_crit + r);

            imgFS(x,y,z) *= amp;
        }

        FFT::inverseFourierTransform(imgFS, img, FFT::Both);

        img.write(fn_out);
    }
	
	return 0;
}
