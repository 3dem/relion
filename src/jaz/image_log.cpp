#include <src/jaz/image_log.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/img_proc/filter_helper.h>


void ImageLog::write(
    Image<Complex> &img, std::string fn,
    bool polar, Centering center,
    double originX, double originY, double originZ,
    double spacingX, double spacingY, double spacingZ)
{
    if (polar)
    {
        Image<RFLOAT> argImg, absImg;
        FilterHelper::getPhase(img, argImg);
        FilterHelper::getAbs(img, absImg);

        write(argImg, fn+"_arg", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
        write(absImg, fn+"_abs", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
    }
    else
    {
        Image<RFLOAT> realImg, imagImg;
        FilterHelper::getReal(img, realImg);
        FilterHelper::getImag(img, imagImg);

        write(realImg, fn+"_re", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
        write(imagImg, fn+"_im", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
    }
}
