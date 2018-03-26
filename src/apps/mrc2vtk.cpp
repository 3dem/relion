#include <src/image.h>
#include <src/jaz/vtk_helper.h>

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "usage: relion_mrc2vtk X\n        X -> X.vtk\n";
        return 666;
    }

    std::string fn(argv[1]);

    Image<RFLOAT> img;
    img.read(fn);

    Image<RFLOAT> imgZ = VtkHelper::allToZ(img);
    VtkHelper::writeVTK(imgZ, fn+".vtk");
}
