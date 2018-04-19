#include <src/jaz/complex_io.h>

void ComplexIO::write(const MultidimArray<Complex>& img, std::string fnBase, std::string fnSuffix)
{
    Image<RFLOAT> temp(img.xdim, img.ydim, img.zdim, img.ndim);

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img)
    {
        DIRECT_NZYX_ELEM(temp.data, l, k, i, j) = DIRECT_NZYX_ELEM(img, l, k, i, j).real;
    }

    temp.write(fnBase + "_real" + fnSuffix);

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img)
    {
        DIRECT_NZYX_ELEM(temp.data, l, k, i, j) = DIRECT_NZYX_ELEM(img, l, k, i, j).imag;
    }

    temp.write(fnBase + "_imag" + fnSuffix);
}

void ComplexIO::read(Image<Complex>& img, std::string fnBase, std::string fnSuffix)
{
    Image<RFLOAT> temp;

    temp.read(fnBase + "_real" + fnSuffix);

    img = Image<Complex>(temp.data.xdim, temp.data.ydim, temp.data.zdim, temp.data.ndim);

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
    {
        DIRECT_NZYX_ELEM(img.data, l, k, i, j).real = DIRECT_NZYX_ELEM(temp.data, l, k, i, j);
    }

    temp.read(fnBase + "_imag" + fnSuffix);

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
    {
        DIRECT_NZYX_ELEM(img.data, l, k, i, j).imag = DIRECT_NZYX_ELEM(temp.data, l, k, i, j);
    }
}
