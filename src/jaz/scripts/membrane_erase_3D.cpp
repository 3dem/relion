
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/membrane/blob_3d.h>
#include <src/jaz/membrane/blob_fit_3d.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/membrane/membrane_fit.h>

#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{
    std::string dir = "/home/zivanovj/oishik/bin8_SIRT_box64_mrc/";

    Catalogue cat("particles/particles.tbl");

    const int pc = cat.particles.size();
    const int num_threads = 16;

	
    #pragma omp parallel for num_threads(num_threads)
    for (int p = 0; p < pc; p++)
    {
        std::cout << p << std::endl;

        DynamoParticle& pp = cat.particles[p];
        std::string index = "000000";

        const int tag = pp.tag;
        if (tag >= 0) index[5] = '0' + tag%10;
        if (tag >= 10) index[4] = '0' + (tag/10)%10;
        if (tag >= 100) index[3] = '0' + (tag/100)%10;
        if (tag >= 1000) index[2] = '0' + (tag/1000)%10;
        if (tag >= 10000) index[1] = '0' + (tag/10000)%10;
        if (tag >= 100000) index[0] = '0' + (tag/100000)%10;

        std::string fn = dir + "particle_" + index + ".mrc";

        try
        {
            BufferedImage<float> box;
            box.read(fn);

            const int w = box.xdim;
            const int h = box.ydim;
            const int d = box.zdim;

            d4Matrix A = pp.getAlignmentMatrixAlias4x4(w,h,d);
            A.invert();

            /*Image<float> box2(w,h,d);

            for (int z = 0; z < d; z++)
            for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
            {
                d4Vector s0(x - w/2, y - h/2, z - d/2, 0.0);
                d4Vector s1 = A * s0 + d4Vector(w/2, h/2, d/2);

                box2(x,y,z) = Interpolation::linearXYZ_clip(box, s1.x, s1.y, s1.z);
            }

            box2.write("debug/"+index+".mrc");*/

            MembraneFit mf(box, 15.0);

            //std::cout << "A = " << A << std::endl;

            d3Vector north(A(0,2), A(1,2), A(2,2));

            //std::cout << "north = " << north << std::endl;

            std::vector<double> initial(9, 0.0);
            initial[3] = north[0];
            initial[6] = north[1];
            initial[8] = north[2];

            BufferedImage<float> reconst0 = mf.expand(initial);
            reconst0.write("debug/"+index+"_reconst0.mrc");

            std::vector<double> opt = NelderMead::optimize(
                initial, mf, 0.001, 1e-8, 1000,
                1, 2, 0.5, 0.5, false);

            BufferedImage<float> reconst = mf.expand(opt);

            reconst.write("membrane_erased/"+index+"_membrane.mrc");
            (box - reconst).write("membrane_erased/"+index+"_rest.mrc");

            std::ofstream os("membrane_erased/"+index+"_coeffs.dat");

            for (int i = 0; i < 9; i++)
            {
                os << opt[i] << ' ';
            }

            os << '\n';

        }
        catch (RelionError e)
        {
        }
    }



    return 0;
}

