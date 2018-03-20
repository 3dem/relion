#include <src/jaz/vtk_helper.h>

using namespace gravis;

Image<RFLOAT> VtkHelper::allToZ(const Image<RFLOAT> &img)
{
    std::cout   << img.data.xdim << "x"
                << img.data.ydim << "x"
                << img.data.zdim << "x"
                << img.data.ndim << "\n";

    if (img.data.ndim == 1) return img;

    Image<RFLOAT> out(img.data.xdim, img.data.ydim, img.data.ndim);

    for (int n = 0; n < img.data.ndim; n++)
    for (int y = 0; y < img.data.ydim; y++)
    for (int x = 0; x < img.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(out(), 0, n, y, x) = DIRECT_NZYX_ELEM(img(), n, 0, y, x);
    }

    return out;
}

void VtkHelper :: writeVTK(Image<double>& img, std::string fn,
                           double originX, double originY, double originZ,
                           double spacingX, double spacingY, double spacingZ,
                           bool binary)
{
    const size_t size = (img.data.xdim * img.data.ydim * img.data.zdim);
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "double";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.data.xdim << " " << img.data.ydim << " " << img.data.zdim << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(img.data.data), sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img.data)
        {
            os << DIRECT_A3D_ELEM(img.data, k, i, j) << "\n";
        }
    }
}

void VtkHelper :: writeVTK(Image<float>& img, std::string fn,
                           double originX, double originY, double originZ,
                           double spacingX, double spacingY, double spacingZ,
                           bool binary)
{
    const size_t size = (img.data.xdim * img.data.ydim * img.data.zdim);
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.data.xdim << " " << img.data.ydim << " " << img.data.zdim << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(img.data.data), sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img.data)
        {
            os << DIRECT_A3D_ELEM(img.data, k, i, j) << "\n";
        }
    }
}

void VtkHelper::writeVTK(Image<Complex> &img, std::string fn, double originX, double originY, double originZ, double spacingX, double spacingY, double spacingZ, bool binary)
{
    const size_t size = (img.data.xdim * img.data.ydim * img.data.zdim);
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    //std::cout << "size: " << size << "\n";
    //std::cout << "sizetype: " << sizetype << "\n";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.data.xdim << " " << img.data.ydim << " " << img.data.zdim << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 2\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(img.data.data), 2*sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img.data)
        {
            os << DIRECT_A3D_ELEM(img.data, k, i, j).real << "\n";
            os << DIRECT_A3D_ELEM(img.data, k, i, j).imag << "\n";
        }
    }
}

void VtkHelper :: writeVTK(MultidimArray<RFLOAT>& img, std::string fn,
                           double originX, double originY, double originZ,
                           double spacingX, double spacingY, double spacingZ,
                           bool binary)
{
    const size_t size = (img.xdim * img.ydim * img.zdim);
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    //std::cout << "size: " << size << "\n";
    //std::cout << "sizetype: " << sizetype << "\n";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.xdim << " " << img.ydim << " " << img.zdim << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(img.data), sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img)
        {
            os << DIRECT_A3D_ELEM(img, k, i, j) << "\n";
        }
    }
}

void VtkHelper :: writeVTK_Complex(MultidimArray<Complex>& img, std::string fn, bool binary)
{
    const size_t size = (img.xdim * img.ydim * img.zdim);
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    //std::cout << "size: " << size << "\n";
    //std::cout << "sizetype: " << sizetype << "\n";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.xdim << " " << img.ydim << " " << img.zdim << "\n";
    os << "SPACING 1 1 1\n";
    os << "ORIGIN 0 0 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 2\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(img.data), 2*sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img)
        {
            os << DIRECT_A3D_ELEM(img, k, i, j).real << "\n";
            os << DIRECT_A3D_ELEM(img, k, i, j).imag << "\n";
        }
    }
}

void VtkHelper :: writeVTK_d3(MultidimArray<gravis::t3Vector<RFLOAT> >& img, std::string fn, bool binary)
{
    const size_t size = (img.xdim * img.ydim * img.zdim);
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    //std::cout << "size: " << size << "\n";
    //std::cout << "sizetype: " << sizetype << "\n";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.xdim << " " << img.ydim << " " << img.zdim << "\n";
    os << "SPACING 1 1 1\n";
    os << "ORIGIN 0 0 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 3\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(img.data), 2*sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img)
        {
            os << DIRECT_A3D_ELEM(img, k, i, j).x << "\n";
            os << DIRECT_A3D_ELEM(img, k, i, j).y << "\n";
            os << DIRECT_A3D_ELEM(img, k, i, j).z << "\n";
        }
    }
}

void VtkHelper :: writeTomoVTK(Image<RFLOAT>& img, std::string fn, bool binary)
{
    const size_t size = (img.data.xdim * img.data.ydim * img.data.ndim);
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    //std::cout << "size: " << size << "\n";
    //std::cout << "sizetype: " << sizetype << "\n";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << img.data.xdim << " " << img.data.ydim << " " << img.data.ndim << "\n";
    os << "SPACING 1 1 1\n";
    os << "ORIGIN 0 0 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(img.data.data), sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
        {
            os << DIRECT_NZYX_ELEM(img.data, l, 0, i, j) << "\n";
        }
    }
}

void VtkHelper :: write(std::vector<Image<RFLOAT> >& stack, std::string fn,
                        double originX, double originY,
                        double spacingX, double spacingY,
                        bool binary)
{
    const size_t size = (stack[0].data.xdim * stack[0].data.ydim * stack.size());
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << stack[0].data.xdim << " " << stack[0].data.ydim << " " << stack.size() << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " 1\n";
    os << "ORIGIN " << originX << " " << originY << " 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(stack[0].data.data), sizeof(RFLOAT)*size);
    }
    else
    {
        for (int ind = 0; ind < stack.size(); ind++)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(stack[ind].data)
            {
                os << DIRECT_A2D_ELEM(stack[ind].data, i, j) << "\n";
            }
        }
    }
}

void VtkHelper :: writeCentered(std::vector<Image<RFLOAT> >& stack, std::string fn,
                        double originX, double originY,
                        double spacingX, double spacingY,
                        bool binary)
{
    const size_t size = (stack[0].data.xdim * stack[0].data.ydim * stack.size());
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << stack[0].data.xdim << " " << stack[0].data.ydim << " " << stack.size() << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " 1\n";
    os << "ORIGIN " << originX << " " << originY << " 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(stack[0].data.data), sizeof(RFLOAT)*size);
    }
    else
    {
        const int w = stack[0].data.xdim;
        const int h = stack[0].data.ydim;

        for (int ind = 0; ind < stack.size(); ind++)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(stack[ind].data)
            {
                int ii = (h + i - h/2)%h;
                int jj = (w + j - w/2)%w;

                os << DIRECT_A2D_ELEM(stack[ind].data, ii, jj) << "\n";
            }
        }
    }
}

void VtkHelper :: write(std::vector<Image<float> >& stack, std::string fn,
                        double originX, double originY,
                        double spacingX, double spacingY,
                        bool binary)
{
    const size_t size = (stack[0].data.xdim * stack[0].data.ydim * stack.size());
    std::ofstream os(fn.c_str(), std::ios::binary);

    std::string sizetype = "float";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << stack[0].data.xdim << " " << stack[0].data.ydim << " " << stack.size() << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " 1\n";
    os << "ORIGIN " << originX << " " << originY << " 0\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(stack[0].data.data), sizeof(RFLOAT)*size);
    }
    else
    {
        for (int ind = 0; ind < stack.size(); ind++)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(stack[ind].data)
            {
                os << DIRECT_A2D_ELEM(stack[ind].data, i, j) << "\n";
            }
        }
    }
}

void VtkHelper :: writeVTK(Volume<RFLOAT>& vol, std::string fn,
                           double originX, double originY, double originZ,
                           double spacingX, double spacingY, double spacingZ,
                           bool binary)
{
    const size_t size = (vol.dimx * vol.dimy * vol.dimz);
    std::ofstream os(fn.c_str());

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << vol.dimx << " " << vol.dimy << " " << vol.dimz << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 1\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(vol.data()), sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_VOXELS(vol)
        {
            os << vol(x,y,z) << "\n";
        }
    }
}

void VtkHelper :: writeVTK(Volume<t3Vector<RFLOAT> >& vol, std::string fn,
                           double originX, double originY, double originZ,
                           double spacingX, double spacingY, double spacingZ,
                           bool binary)
{
    const size_t size = (vol.dimx * vol.dimy * vol.dimz);
    std::ofstream os(fn.c_str());

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << vol.dimx << " " << vol.dimy << " " << vol.dimz << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 3\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(vol.data()), sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_VOXELS(vol)
        {
            os << vol(x,y,z).x << " " << vol(x,y,z).y << " " << vol(x,y,z).z << "\n";
        }
    }
}

void VtkHelper :: readVTK(std::string fn, Volume<RFLOAT>& vol, d3Vector& origin, d3Vector& spacing)
{
    std::ifstream file(fn.c_str());

    char text[4096];

    if (!file.is_open())
    {
        REPORT_ERROR("failed to open " + fn + '\n');
    }

    file.getline(text, 4096);

    if (std::string(text) != "# vtk DataFile Version 2.0")
    {
        REPORT_ERROR("Unsupported VTK format: " + std::string(text) + '\n');
    }

    file.getline(text, 4096);
    file.getline(text, 4096);

    if (std::string(text) != "ASCII")
    {
        REPORT_ERROR("Only ASCII VTK files are supported.\n");
    }
    /*
    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << vol.dimx << " " << vol.dimy << " " << vol.dimz << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 6\n";
    os << "LOOKUP_TABLE default\n";
    */

    file.getline(text, 4096);

    std::string first, dummy;
    int dimx, dimy, dimz, dims;
    size_t size;

    for (int i = 0; i < 6; i++)
    {
        file.getline(text, 4096);
        std::stringstream linestream(text);

        linestream >> first;

        if (first == "DIMENSIONS")
        {
            linestream >> dimx;
            linestream >> dimy;
            linestream >> dimz;

        }
        else if (first == "SPACING")
        {
            linestream >> spacing.x;
            linestream >> spacing.y;
            linestream >> spacing.z;
        }
        else if (first == "ORIGIN")
        {
            linestream >> origin.x;
            linestream >> origin.y;
            linestream >> origin.z;
        }
        else if (first == "POINT_DATA")
        {
            linestream >> size;
        }
        else if (first == "SCALARS")
        {
            linestream >> dummy;
            linestream >> dummy;
            linestream >> dims;

            if (dims != 1)
            {
                std::stringstream sts;
                std::string st;
                sts << fn << " is not a scalar volume (voxeldims = " << dims << ")\n";
                sts >> st;

                REPORT_ERROR(st);
            }
        }
        else if (first == "LOOKUP_TABLE")
        {
            linestream >> dummy;
        }
    }

    if (size != ((size_t)dimx)*((size_t)dimy)*((size_t)dimz))
    {
        std::cout << "Bad size info in " << fn << ": " << size << " vs. " << (dimx*dimy*dimz) << "\n";
        std::exit(666);

        std::stringstream sts;
        std::string st;
        sts << "Bad size info in " << fn << ": " << size << " vs. " << (dimx*dimy*dimz) << "\n";
        sts >> st;

        REPORT_ERROR(st);
    }

    vol.resize(dimx, dimy, dimz);

    for (size_t i = 0; i < size; i++)
    {
        file >> vol.voxels[i];
    }
}

void VtkHelper :: writeVTK(Volume<Tensor3x3<RFLOAT> >& vol, std::string fn,
                           double originX, double originY, double originZ,
                           double spacingX, double spacingY, double spacingZ,
                           bool binary)
{
    const size_t size = (vol.dimx * vol.dimy * vol.dimz);
    std::ofstream os(fn.c_str());

    std::string sizetype = "float";
    if (sizeof(RFLOAT) > 4) sizetype = "double";

    std::cout << "writing: " << fn << "\n";

    os << "# vtk DataFile Version 2.0\n";
    os << "Volume example\n";

    if (binary)
    {
        os << "BINARY\n";
    }
    else
    {
        os << "ASCII\n";
    }

    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << vol.dimx << " " << vol.dimy << " " << vol.dimz << "\n";
    os << "SPACING " << spacingX << " " << spacingY << " " << spacingZ << "\n";
    os << "ORIGIN " << originX << " " << originY << " " << originZ << "\n";
    os << "POINT_DATA " << size << "\n";
    os << "SCALARS volume_scalars " << sizetype << " 6\n";
    os << "LOOKUP_TABLE default\n";

    if (binary)
    {
        os.write((char*)(vol.data()), sizeof(RFLOAT)*size);
    }
    else
    {
        FOR_ALL_VOXELS(vol)
        {
            os << vol(x,y,z).xx << " " << vol(x,y,z).yy << " " << vol(x,y,z).zz << " " << vol(x,y,z).xy << " " << vol(x,y,z).yz << " " << vol(x,y,z).xz << "\n";
        }
    }
}

