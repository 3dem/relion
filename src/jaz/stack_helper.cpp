#include <src/jaz/stack_helper.h>
#include <src/jaz/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/Fourier_helper.h>
#include <src/jaz/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/vtk_helper.h>
#include <src/fftw.h>

#include <omp.h>

using namespace gravis;

std::vector<MetaDataTable> StackHelper::splitByStack(const MetaDataTable* mdt)
{
    std::vector<MetaDataTable> out(0);

    if (!mdt->labelExists(EMDL_IMAGE_NAME))
    {
        REPORT_ERROR("StackHelper::splitByStack: "+EMDL::label2Str(EMDL_IMAGE_NAME)+" missing in meta_data_table.\n");
    }

    std::string testString;
    mdt->getValue(EMDL_IMAGE_NAME, testString, 0);

    if (testString.find("@") < 0)
    {
        REPORT_ERROR("StackHelper::splitByStack: "+EMDL::label2Str(EMDL_IMAGE_NAME)+" does not contain an '@'.\n");
    }

    MetaDataTable md2(*mdt);
    md2.newSort(EMDL_IMAGE_NAME, false, true);

    const long lc = md2.numberOfObjects();
    std::string lastName = "", curName, curFullName;
    long curInd = -1;

    for (int i = 0; i < lc; i++)
    {
        md2.getValue(EMDL_IMAGE_NAME, curFullName, i);

        curName = curFullName.substr(curFullName.find("@")+1);

        if (curName != lastName)
        {
            lastName = curName;
            curInd++;
            out.push_back(MetaDataTable());
        }

        out[curInd].addObject(md2.getObject(i));
    }

    for (int i = 0; i <= curInd; i++)
    {
        out[i].newSort(EMDL_IMAGE_NAME, false, false, true);
    }

    return out;
}

std::vector<Image<RFLOAT> > StackHelper::loadStack(const MetaDataTable* mdt, std::string path, int threads)
{
    std::vector<Image<RFLOAT>> out(mdt->numberOfObjects());
    const long ic = mdt->numberOfObjects();

    std::string name, fullName;
    mdt->getValue(EMDL_IMAGE_NAME, fullName, 0);
    name = fullName.substr(fullName.find("@")+1);

    if (path != "")
    {
        name = path + "/" + name.substr(name.find_last_of("/")+1);
    }

    Image<RFLOAT> in;
    in.read(name);

    const int w = in.data.xdim;
    const int h = in.data.ydim;

    #pragma omp parallel for num_threads(threads)
    for (long i = 0; i < ic; i++)
    {
        std::string sliceName;
        mdt->getValue(EMDL_IMAGE_NAME, sliceName, i);
        std::string indStr = sliceName.substr(0,sliceName.find("@"));

        std::istringstream str(indStr);
        int j;
        str >> j;

        out[i] = Image<RFLOAT>(w,h);
        SliceHelper::extractStackSlice(in, out[i], j-1);
    }

    return out;
}

std::vector<Image<Complex> > StackHelper::loadStackFS(const MetaDataTable *mdt, std::string path,
                                                     int threads, std::vector<FourierTransformer> *fts)
{
    std::vector<Image<Complex> > out(mdt->numberOfObjects());
    const long ic = mdt->numberOfObjects();

    std::string name, fullName;
    mdt->getValue(EMDL_IMAGE_NAME, fullName, 0);
    name = fullName.substr(fullName.find("@")+1);

    if (path != "")
    {
        name = path + "/" + name.substr(name.find_last_of("/")+1);
    }

    Image<RFLOAT> in;
    in.read(name);

    const int w = in.data.xdim;
    const int h = in.data.ydim;

    #pragma omp parallel for num_threads(threads)
    for (long i = 0; i < ic; i++)
    {
        int threadnum = omp_get_thread_num();

        std::string sliceName;
        mdt->getValue(EMDL_IMAGE_NAME, sliceName, i);
        std::string indStr = sliceName.substr(0,sliceName.find("@"));

        std::istringstream str(indStr);
        int j;
        str >> j;

        Image<RFLOAT> aux(w,h);
        SliceHelper::extractStackSlice(in, aux, j-1);
        (*fts)[threadnum].FourierTransform(aux(), out[i]());
    }

    return out;
}

void StackHelper::saveStack(std::vector<Image<RFLOAT> > &stack, std::string fn)
{
    const int w = stack[0].data.xdim;
    const int h = stack[0].data.ydim;
    const int c = stack.size();

    Image<RFLOAT> img(w,h,1,c);

    for (int i = 0; i < c; i++)
    {
        SliceHelper::insertStackSlice(stack[i], img, i);
    }

    img.write(fn);
}

std::vector<std::vector<Image<RFLOAT> > > StackHelper::loadMovieStack(const MetaDataTable* mdt, std::string moviePath)
{
    std::vector<std::vector<Image<RFLOAT> > > out(mdt->numberOfObjects());
    const long pc = mdt->numberOfObjects();

    std::string name, fullName, movieName;
    mdt->getValue(EMDL_IMAGE_NAME, fullName, 0);
    mdt->getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
    name = fullName.substr(fullName.find("@")+1);

    std::string finName;

    if (moviePath == "")
    {
        finName = name;
    }
    else
    {
        finName = moviePath + "/" + movieName.substr(movieName.find_last_of("/")+1);
    }

    std::cout << "loading real: " << finName << "\n";

    Image<RFLOAT> in;
    in.read(finName);

    std::cout << "size = " << in.data.xdim << "x" << in.data.ydim << "x" << in.data.zdim << "x" << in.data.ndim << "\n";
    std::cout << "pc = " << pc << "\n";

    const int fc = in.data.ndim / pc;

    const int w = in.data.xdim;
    const int h = in.data.ydim;

    for (long p = 0; p < pc; p++)
    {
        out[p] = std::vector<Image<RFLOAT> >(fc);

        for (long f = 0; f < fc; f++)
        {
            out[p][f] = Image<RFLOAT>(w,h);
            SliceHelper::extractStackSlice(in, out[p][f], f*pc + p);
        }
    }

    return out;
}

std::vector<std::vector<Image<Complex> > > StackHelper::loadMovieStackFS(const MetaDataTable* mdt,
                                                                         std::string moviePath,
                                                                         bool center, int threads,
                                                                         std::vector<FourierTransformer> *fts)
{
    std::vector<std::vector<Image<Complex> > > out(mdt->numberOfObjects());
    const long pc = mdt->numberOfObjects();

    std::string name, fullName, movieName;
    mdt->getValue(EMDL_IMAGE_NAME, fullName, 0);
    mdt->getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
    name = fullName.substr(fullName.find("@")+1);

    std::string finName;

    if (moviePath == "")
    {
        finName = name;
    }
    else
    {
        finName = moviePath + "/" + movieName.substr(movieName.find_last_of("/")+1);
    }

    std::cout << "loading " << finName << "\n";

    Image<RFLOAT> in;
    in.read(finName);

    std::cout << "size = " << in.data.xdim << "x" << in.data.ydim << "x" << in.data.zdim << "x" << in.data.ndim << "\n";
    std::cout << "pc = " << pc << "\n";


    const int fc = in.data.ndim / pc;

    const int w = in.data.xdim;
    const int h = in.data.ydim;

    for (long p = 0; p < pc; p++)
    {
        out[p] = std::vector<Image<Complex> >(fc);

        if (threads > 1)
        {
            #pragma omp parallel for num_threads(threads)
            for (long f = 0; f < fc; f++)
            {
                int threadnum = omp_get_thread_num();

                Image<RFLOAT> aux(w,h);

                SliceHelper::extractStackSlice(in, aux, f*pc + p);
                if (center) CenterFFT(aux(), true);
                (*fts)[threadnum].FourierTransform(aux(), out[p][f]());
            }
        }
        else
        {
            Image<RFLOAT> aux(w,h);
            FourierTransformer ft;

            for (long f = 0; f < fc; f++)
            {
                SliceHelper::extractStackSlice(in, aux, f*pc + p);
                if (center) CenterFFT(aux(), true);
                ft.FourierTransform(aux(), out[p][f]());
            }
        }
    }

    return out;
}

Image<Complex> StackHelper::projectView(Projector* projector, const MetaDataTable* mdt, int index)
{
    const int s = projector->ori_size;
    const int sh = s/2 + 1;

    double xoff, yoff;

    mdt->getValue(EMDL_ORIENT_ORIGIN_X, xoff, index);
    mdt->getValue(EMDL_ORIENT_ORIGIN_Y, yoff, index);

    double rot, tilt, psi;

    Matrix2D<RFLOAT> A3D;
    mdt->getValue(EMDL_ORIENT_ROT, rot, index);
    mdt->getValue(EMDL_ORIENT_TILT, tilt, index);
    mdt->getValue(EMDL_ORIENT_PSI, psi, index);
    Euler_angles2matrix(rot, tilt, psi, A3D);

    Image<Complex> out(sh, s);
    out.data.initZeros();

    projector->get2DFourierTransform(out.data, A3D, false);
    shiftImageInFourierTransform(out(), out(), s, s/2 - xoff, s/2 - yoff);

    return out;
}

std::vector<Image<Complex> > StackHelper::projectStack(Projector* projector, const MetaDataTable* mdt)
{
    std::vector<Image<Complex> > out(mdt->numberOfObjects());
    const long ic = mdt->numberOfObjects();

    const int s = projector->ori_size;
    const int sh = s/2 + 1;

    for (long i = 0; i < ic; i++)
    {
        double xoff, yoff;

        mdt->getValue(EMDL_ORIENT_ORIGIN_X, xoff, i);
        mdt->getValue(EMDL_ORIENT_ORIGIN_Y, yoff, i);

        double rot, tilt, psi;

        Matrix2D<RFLOAT> A3D;
        mdt->getValue(EMDL_ORIENT_ROT, rot, i);
        mdt->getValue(EMDL_ORIENT_TILT, tilt, i);
        mdt->getValue(EMDL_ORIENT_PSI, psi, i);
        Euler_angles2matrix(rot, tilt, psi, A3D);

        out[i] = Image<Complex>(sh, s);
        out[i].data.initZeros();
        projector->get2DFourierTransform(out[i].data, A3D, false);
        shiftImageInFourierTransform(out[i](), out[i](), s, s/2 - xoff, s/2 - yoff);
    }

    return out;
}

std::vector<Image<Complex> > StackHelper::projectStackPar(
        Projector* projector, const MetaDataTable* mdt, int numThreads)
{
    std::vector<Image<Complex> > out(mdt->numberOfObjects());
    const long ic = mdt->numberOfObjects();

    const int s = projector->ori_size;
    const int sh = s/2 + 1;

    #pragma omp parallel for num_threads(numThreads)
    for (long i = 0; i < ic; i++)
    {
        double xoff, yoff;

        mdt->getValue(EMDL_ORIENT_ORIGIN_X, xoff, i);
        mdt->getValue(EMDL_ORIENT_ORIGIN_Y, yoff, i);

        double rot, tilt, psi;

        Matrix2D<RFLOAT> A3D;
        mdt->getValue(EMDL_ORIENT_ROT, rot, i);
        mdt->getValue(EMDL_ORIENT_TILT, tilt, i);
        mdt->getValue(EMDL_ORIENT_PSI, psi, i);
        Euler_angles2matrix(rot, tilt, psi, A3D);

        out[i] = Image<Complex>(sh, s);
        out[i].data.initZeros();
        projector->get2DFourierTransform(out[i].data, A3D, false);
        shiftImageInFourierTransform(out[i](), out[i](), s, s/2 - xoff, s/2 - yoff);
    }

    return out;
}

std::vector<Image<Complex> > StackHelper::FourierTransform(std::vector<Image<RFLOAT> >& stack)
{
    std::vector<Image<Complex> > out(stack.size());
    const long ic = stack.size();

    for (long i = 0; i < ic; i++)
    {
        FourierTransformer ft;
        ft.FourierTransform(stack[i].data, out[i].data);
    }

    return out;
}

std::vector<Image<RFLOAT> > StackHelper::inverseFourierTransform(std::vector<Image<Complex> >& stack)
{
    std::vector<Image<RFLOAT> > out(stack.size());
    const long ic = stack.size();

    std::cout << ic << " images in stack.\n";

    const int h = stack[0].data.ydim;
    const int ww = stack[0].data.xdim;
    const int w = 2*(ww - 1);

    for (long i = 0; i < ic; i++)
    {
        out[i] = Image<RFLOAT>(w,h);

        FourierTransformer ft;
        ft.inverseFourierTransform(stack[i].data, out[i].data);
    }

    return out;
}

std::vector<Image<Complex> > StackHelper::applyBeamTilt(std::vector<Image<Complex> > &stack,
                                RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                                RFLOAT tilt_x, RFLOAT tilt_y)
{
    std::vector<Image<Complex> > out(stack.size());
    const long ic = stack.size();

    for (long i = 0; i < ic; i++)
    {
        out[i] = stack[i];
        selfApplyBeamTilt(out[i].data, tilt_x, tilt_y, lambda, Cs, angpix, stack[i].data.ydim);
    }

    return out;

}

std::vector<Image<Complex> > StackHelper::applyBeamTiltPar(std::vector<Image<Complex> > &stack,
                                RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                                RFLOAT tilt_x, RFLOAT tilt_y, int numThreads)
{
    std::vector<Image<Complex> > out(stack.size());
    const long ic = stack.size();

    #pragma omp parallel for num_threads(numThreads)
    for (long i = 0; i < ic; i++)
    {
        out[i] = stack[i];
        selfApplyBeamTilt(out[i].data, tilt_x, tilt_y, lambda, Cs, angpix, stack[i].data.ydim);
    }

    return out;

}

std::vector<double> StackHelper::powerSpectrum(const std::vector<std::vector<Image<Complex>>> &stack)
{
    const int ic = stack.size();
    const int fc = stack[0].size();
    const int w = stack[0][0].data.xdim;
    const int h = stack[0][0].data.ydim;

    std::vector<double> out(w, 0.0), wgh(w, 0.0);

    for (int i = 0; i < ic; i++)
    for (int f = 0; f < fc; f++)
    {
        for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            const Complex z = DIRECT_A2D_ELEM(stack[i][f].data, y, x);

            const double yy = y < w? y : y - h;
            const double xx = x;

            const int r = (int) sqrt(xx*xx + yy*yy);

            if (r >= w) continue;

            out[r] += z.norm();
            wgh[r] += 1.0;
        }
    }

    for (int x = 0; x < w; x++)
    {
        if (wgh[x] > 0.0)
        {
            out[x] /= wgh[x];
        }
    }

    return out;
}

std::vector<double> StackHelper::varSpectrum(const std::vector<std::vector<Image<Complex>>> &stack)
{
    const int ic = stack.size();
    const int fc = stack[0].size();
    const int w = stack[0][0].data.xdim;
    const int h = stack[0][0].data.ydim;

    std::vector<double> out(w, 0.0), wgh(w, 0.0);
    std::vector<Complex> mean(w, 0.0);

    for (int i = 0; i < ic; i++)
    for (int f = 0; f < fc; f++)
    {
        for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            const Complex z = DIRECT_A2D_ELEM(stack[i][f].data, y, x);

            const double yy = y < w? y : y - h;
            const double xx = x;

            const int r = (int) sqrt(xx*xx + yy*yy);

            if (r >= w) continue;

            mean[r] += z;
            wgh[r] += 1.0;
        }
    }

    for (int x = 0; x < w; x++)
    {
        if (wgh[x] > 0.0)
        {
            mean[x] /= wgh[x];
        }
    }

    for (int i = 0; i < ic; i++)
    for (int f = 0; f < fc; f++)
    {
        for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            const Complex z = DIRECT_A2D_ELEM(stack[i][f].data, y, x);

            const double yy = y < w? y : y - h;
            const double xx = x;

            const int r = (int) sqrt(xx*xx + yy*yy);

            if (r >= w) continue;

            out[r] += (z - mean[r]).norm();
        }
    }

    for (int x = 0; x < w; x++)
    {
        if (wgh[x] > 1.0)
        {
            out[x] /= (wgh[x] - 1.0);
        }
    }

    return out;
}

std::vector<double> StackHelper::powerSpectrum(
        const std::vector<std::vector<Image<Complex>>>& obs,
        const std::vector<Image<Complex> >& signal)
{
    const int ic = obs.size();
    const int fc = obs[0].size();
    const int w = obs[0][0].data.xdim;
    const int h = obs[0][0].data.ydim;

    std::vector<double> out(w, 0.0), wgh(w, 0.0);

    for (int i = 0; i < ic; i++)
    for (int f = 0; f < fc; f++)
    {
        for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            const Complex z = DIRECT_A2D_ELEM(obs[i][f].data, y, x) - DIRECT_A2D_ELEM(signal[i].data, y, x);

            const double yy = y < w? y : y - h;
            const double xx = x;

            const int r = (int) sqrt(xx*xx + yy*yy);

            if (r >= w) continue;

            out[r] += z.norm();
            wgh[r] += 1.0;
        }
    }

    for (int x = 0; x < w; x++)
    {
        if (wgh[x] > 0.0)
        {
            out[x] /= wgh[x];
        }
    }

    return out;
}
