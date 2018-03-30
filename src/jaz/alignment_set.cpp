#include <src/jaz/alignment_set.h>
#include <src/jaz/vtk_helper.h>

#include <omp.h>

using namespace gravis;

AlignmentSet::AlignmentSet()
:   mc(0), fc(0), s(0), sh(0), k0(0), k1(0)
{
}

AlignmentSet::AlignmentSet(
        const std::vector<MetaDataTable> &mdts,
        int fc, int s, int k0, int k1)
:   mc(mdts.size()),
    fc(fc),
    s(s),
    sh(s/2+1),
    k0(k0),
    k1(k1)
{
    CCs.resize(mc);
    obs.resize(mc);
    pred.resize(mc);

    positions.resize(mc);
    initialTracks.resize(mc);
    globComp.resize(mc);

    for (int m = 0; m < mc; m++)
    {
        const int pc = mdts[m].numberOfObjects();
        CCs[m] = std::vector<std::vector<Image<double>>>(pc, std::vector<Image<double>>(fc));
        obs[m] = std::vector<std::vector<std::vector<d2Vector>>>(pc, std::vector<std::vector<d2Vector>>(fc));
        pred[m].resize(pc);
    }

    accCoords.reserve(sh*s);

    int num = 0;

    for (int y = 0; y < s; y++)
    for (int x = 0; x < sh; x++)
    {
        const double xx = x;
        const double yy = y < sh? y : y - s;

        int r = ROUND(sqrt(xx*xx + yy*yy));

        if (r >= k0 && r < k1)
        {
            accCoords.push_back(t2Vector<int>(x,y));
            num++;
        }
    }

    accPix = num;
}

std::vector<d2Vector> AlignmentSet::accelerate(const Image<Complex> &img)
{
    std::vector<d2Vector> out(accPix);

    for (int i = 0; i < accPix; i++)
    {
        t2Vector<int> c = accCoords[i];

        Complex z = img(c.y, c.x);

        out[i] = d2Vector(z.real, z.imag);
    }

    return out;
}

d3Vector AlignmentSet::updateTsc(
    const std::vector<std::vector<d2Vector>>& tracks,
    int mg, int threads)
{
    const int pad = 512;
    std::vector<d3Vector> outT(pad*threads, d3Vector(0.0, 0.0, 0.0));

    const int pc = tracks.size();

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        int t = omp_get_thread_num();

        const d2Vector shift = tracks[p][f] / s;

        for (int i = 0; i < accPix; i++)
        {
            t2Vector<int> acc = accCoords[i];

            double x = acc.x;
            double y = acc.y < sh? acc.y : acc.y - s;

            const double dotp = 2 * PI * (x * shift.x + y * shift.y);

            double a, b;

            SINCOS(dotp, &b, &a);

            const d2Vector z_obs_f2 = obs[mg][p][f][i];
            const double c = z_obs_f2.x;
            const double d = z_obs_f2.y;

            const double ac = a * c;
            const double bd = b * d;

            const double ab_cd = (a + b) * (c + d);

            const Complex z_obs = Complex(ac - bd, ab_cd - ac - bd);

            const d2Vector z_pred_f2 = pred[mg][p][i];

            const Complex z_pred = Complex(z_pred_f2.x, z_pred_f2.y);

            outT[pad*t][0] += z_pred.real * z_obs.real + z_pred.imag * z_obs.imag;
            outT[pad*t][1] += z_obs.norm();
            outT[pad*t][2] += z_pred.norm();
        }
    }

    d3Vector out(0.0, 0.0, 0.0);

    for (int t = 0; t < threads; t++)
    {
        out += outT[pad*t];
    }

    return out;
}
