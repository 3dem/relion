#include <src/jaz/stack_helper.h>
#include <src/jaz/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/Fourier_helper.h>
#include <src/jaz/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/vtk_helper.h>
#include <src/fftw.h>
#include <src/micrograph_model.h>
#include <src/jaz/resampling_helper.h>

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

std::vector<std::vector<Image<Complex>>> StackHelper::extractMovieStackFS(
    const MetaDataTable *mdt,
    std::string metaPath, std::string moviePath,
    int outBin, int coordsBin, int movieBin, int squareSize,
    std::vector<FourierTransformer>& fts,
    bool useGain, BinningType binningType,
    bool loadData, bool verbose)
{
    std::vector<std::vector<Image<Complex>>> out(mdt->numberOfObjects());
    const long pc = mdt->numberOfObjects();

    std::string name0;
    mdt->getValue(EMDL_MICROGRAPH_NAME, name0, 0);

    std::string ending = name0.substr(name0.find_last_of(".")+1);

    if (verbose)
    {
        std::cout << "orig. name: " << name0 << "\n";
        std::cout << "ending: " << ending << "\n";
    }

    Image<RFLOAT> mgStack;

    if (ending == "star")
    {
        std::string metaNameAct;

        if (metaPath == "")
        {
            metaNameAct = name0;
        }
        else
        {
            metaNameAct = metaPath + "/" + name0.substr(name0.find_last_of("/")+1);
        }

        if (verbose)
        {
            std::cout << "loading: " << metaNameAct << "\n";
        }

        MetaDataTable metaMdt;
        metaMdt.read(metaNameAct, "general");

        std::string mgName0;

        if (!metaMdt.containsLabel(EMDL_MICROGRAPH_MOVIE_NAME))
        {
            REPORT_ERROR("StackHelper::extractMovieStackFS: " + metaNameAct
                         + " does not contain a rlnMicrographMovieName field.");
        }

        metaMdt.getValue(EMDL_MICROGRAPH_MOVIE_NAME, mgName0, 0);

        std::string mgNameAct, gainNameAct;

        if (moviePath == "")
        {
            mgNameAct = mgName0;
        }
        else
        {
            mgNameAct = moviePath + "/" + mgName0.substr(mgName0.find_last_of("/")+1);
        }

        if (verbose)
        {
            std::cout << "loading: " << mgNameAct << "\n";
        }

        mgStack.read(mgNameAct, loadData);

        if (useGain)
        {
            if (!metaMdt.containsLabel(EMDL_MICROGRAPH_GAIN_NAME))
            {
                std::cerr << "warning: " << metaNameAct
                          << " does not contain a rlnMicrographGainName field.\n";
            }
            else
            {
                std::string gainName0;
                metaMdt.getValue(EMDL_MICROGRAPH_GAIN_NAME, gainName0, 0);

                if (moviePath == "")
                {
                    gainNameAct = gainName0;
                }
                else
                {
                    gainNameAct = moviePath + "/" + gainName0.substr(gainName0.find_last_of("/")+1);
                }

                Image<RFLOAT> gainRef;
                gainRef.read(gainNameAct, loadData);

                if (loadData)
                {
                    if (gainRef.data.xdim != mgStack.data.xdim
                        || gainRef.data.ydim != mgStack.data.ydim)
                    {
                        std::stringstream stsm;
                        stsm << mgStack.data.xdim << "x" << mgStack.data.ydim;
                        std::stringstream stsg;
                        stsg << gainRef.data.xdim << "x" << gainRef.data.ydim;

                        REPORT_ERROR("StackHelper::extractMovieStackFS: gain reference ("
                                     + gainNameAct + ") and movie stack (" + mgNameAct
                                     + ") are of unequal size ("+stsg.str()+" vs. "+stsm.str()+").");
                    }

                    for (int n = 0; n < mgStack().ndim; n++)
                    for (int z = 0; z < mgStack().zdim; z++)
                    for (int y = 0; y < mgStack().ydim; y++)
                    for (int x = 0; x < mgStack().xdim; x++)
                    {
                        DIRECT_NZYX_ELEM(mgStack(), n, z, y, x)
                            /= DIRECT_NZYX_ELEM(gainRef(), 0, 0, y, x);
                    }
                }
            }
        }
    }
    else
    {
        if (useGain)
        {
            std::cerr << "Warning: unable to load gain reference - rlnMicrographName (" << name0 << ") is not a .star file.\n";
        }

        std::string mgNameAct;

        if (moviePath == "")
        {
            mgNameAct = name0;
        }
        else
        {
            mgNameAct = moviePath + "/" + name0.substr(name0.find_last_of("/")+1);
        }

        if (verbose)
        {
            std::cout << "loading: " << mgNameAct << "\n";
        }

        mgStack.read(mgNameAct, loadData);
    }

    if (verbose)
    {
        std::cout << "size: "
            << mgStack().xdim << "x"
            << mgStack().ydim << "x"
            << mgStack().zdim << "x"
            << mgStack().ndim << "\n";
    }

    const bool dataInZ = mgStack.data.zdim > 1;

    const int w0 = mgStack.data.xdim;
    const int h0 = mgStack.data.ydim;
    const int fc = dataInZ? mgStack.data.zdim : mgStack.data.ndim;

    const int nsc = dataInZ? 0 : 1;
    const int zsc = dataInZ? 1 : 0;

    if (verbose)
    {
        if (dataInZ) std::cout << "data in Z\n";
        else std::cout << "data in N\n";

        std::cout << "frame count = " << fc << "\n";
    }

    const int sqMg = squareSize * outBin / movieBin;

    for (long p = 0; p < pc; p++)
    {
        out[p] = std::vector<Image<Complex>>(fc);

        if (!loadData) continue;

        double xp0, yp0;

        mdt->getValue(EMDL_IMAGE_COORD_X, xp0, p);
        mdt->getValue(EMDL_IMAGE_COORD_Y, yp0, p);

        const double xpm = xp0 * coordsBin / (double)movieBin;
        const double ypm = yp0 * coordsBin / (double)movieBin;

        const int threads = fts.size();

        #pragma omp parallel for num_threads(threads)
        for (long f = 0; f < fc; f++)
        {
            int threadnum = omp_get_thread_num();

            Image<RFLOAT> aux0(sqMg, sqMg), aux1(squareSize, squareSize);
            Image<Complex> aux2;

            // @TODO: read whole-image shifts from micrograph class

            const int x0 = (int)xpm;
            const int y0 = (int)ypm;

            for (long int y = 0; y < sqMg; y++)
            for (long int x = 0; x < sqMg; x++)
            {
                int xx = x0 + x;
                int yy = y0 + y;

                if (xx < 0) xx = 0;
                else if (xx >= w0) xx = w0 - 1;

                if (yy < 0) yy = 0;
                else if (yy >= h0) yy = h0 - 1;

                DIRECT_NZYX_ELEM(aux0.data, 0, 0, y, x)
                    = DIRECT_NZYX_ELEM(mgStack.data, nsc*f, zsc*f, yy, xx);
            }

            if (outBin == movieBin)
            {
                fts[threadnum].FourierTransform(aux0(), out[p][f]());
            }
            else if (binningType != FourierCrop)
            {
                if (binningType == BoxBin)
                {
                    ResamplingHelper::downsampleBox2D(aux0, outBin/(double)movieBin, aux1);

                    if (p == 0 && f < 2)
                    {
                        std::stringstream sts;
                        sts << p << "_" << f;
                        VtkHelper::writeVTK(aux0, "debug/aux0_"+sts.str()+".vtk");
                        VtkHelper::writeVTK(aux1, "debug/aux1_"+sts.str()+".vtk");
                    }
                }
                else if (binningType == GaussBin)
                {
                    ResamplingHelper::downsampleGauss2D(aux0, outBin/(double)movieBin, aux1);
                }

                fts[threadnum].FourierTransform(aux1(), out[p][f]());
            }
            else
            {
                fts[threadnum].FourierTransform(aux0(), aux2());

                out[p][f] = FilterHelper::cropCorner2D(aux2, squareSize, squareSize);
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
