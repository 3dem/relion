
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

    Catalogue cat0("mask8_SIRT/results/ite_0006/starting_values/starting_table_ref_001_ite_0006.tbl");
    Catalogue cat1("mask9_SIRT/results/ite_0006/starting_values/starting_table_ref_001_ite_0006.tbl");

    const int pc = cat0.particles.size();

    std::cout << pc << " particles found." << std::endl;

    const int w = 64;
    const int h = 64;
    const int d = 64;

    std::ofstream os0("membrane_dist_01.dat");
    std::ofstream os1("membrane_dist_10.dat");
    std::ofstream os("membrane_dist.dat");

    const double d0 = -1.0;
    const double d1 = 1.0;

    Catalogue catOut0_8, catOut1_8, catOut2_8;
    Catalogue catOut0_9, catOut1_9, catOut2_9;

    std::string fn0_8 = "particles_8_by_membrane_dist_0.tbl";
    std::string fn1_8 = "particles_8_by_membrane_dist_1.tbl";
    std::string fn2_8 = "particles_8_by_membrane_dist_2.tbl";

    std::string fn0_9 = "particles_9_by_membrane_dist_0.tbl";
    std::string fn1_9 = "particles_9_by_membrane_dist_1.tbl";
    std::string fn2_9 = "particles_9_by_membrane_dist_2.tbl";


    for (int p = 0; p < pc; p++)
    {
        DynamoParticle pp0 = cat0.particles[p];
        DynamoParticle pp1 = cat1.particles[p];

        d4Matrix A0 = pp0.getAlignmentMatrixAlias4x4(w,h,d);
        d4Matrix A1 = pp1.getAlignmentMatrixAlias4x4(w,h,d);

        A0.invert();
        A1.invert();

        d3Vector z0(A0(0,2), A0(1,2), A0(2,2));
        d3Vector z1(A1(0,2), A1(1,2), A1(2,2));

        d3Vector z0_box(pp0.dx, pp0.dy, pp0.dz);
        d3Vector z1_box(pp1.dx, pp1.dy, pp1.dz);
        d3Vector dz_box = z1_box - z0_box;

        double dz_ref_0 = z0.dot(dz_box);
        double dz_ref_1 = z1.dot(dz_box);
        double dz_ref = (dz_ref_0 + dz_ref_1)/2.0;

        os0 << dz_ref_0 << '\n';
        os1 << dz_ref_1 << '\n';
        os << dz_ref << '\n';

        if (dz_ref < d0)
        {
            catOut0_8.particles.push_back(pp0);
            catOut0_9.particles.push_back(pp1);
        }
        else if (dz_ref < d1)
        {
            catOut1_8.particles.push_back(pp0);
            catOut1_9.particles.push_back(pp1);
        }
        else
        {
            catOut2_8.particles.push_back(pp0);
            catOut2_9.particles.push_back(pp1);
        }
    }

    catOut0_8.write(fn0_8);
    catOut1_8.write(fn1_8);
    catOut2_8.write(fn2_8);

    catOut0_9.write(fn0_9);
    catOut1_9.write(fn1_9);
    catOut2_9.write(fn2_9);

    return 0;
}

