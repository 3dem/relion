
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
    std::string refFn = "/home/zivanovj/oishik/mask4_SIRT/results/ite_0003/averages/average_ref_001_ite_0003.mrc";

    Catalogue cat("/home/zivanovj/oishik/mask4_SIRT/results/ite_0004/starting_values/starting_table_ref_001_ite_0004.tbl");

    const int pc = cat.particles.size();
    const int num_threads = 16;

    BufferedImage<float> ref;
    ref.read(refFn);

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

            BufferedImage<float> box2(w,h,d);

            for (int z = 0; z < d; z++)
            for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
            {
                d4Vector s0(x - w/2, y - h/2, z - d/2, 0.0);
                d4Vector s1 = A * s0 + d4Vector(w/2, h/2, d/2);

                box2(x,y,z) = box(x,y,z) - Interpolation::linearXYZ_clip(ref, s1.x, s1.y, s1.z);
            }

            box2.write("avg_erased/particle_"+index+".mrc");

            std::cout << "written " << ("avg_erased/particle_"+index+".mrc") << std::endl;

        }
        catch (RelionError e)
        {
            std::cout << "unable to read " << fn << std::endl;
        }
    }

    return 0;
}

