/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef ALIGNMENT_SET_H
#define ALIGNMENT_SET_H

#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t3Vector.h>
#include <vector>
#include <omp.h>


template<class T>
class AlignmentSet
{
    public:

        AlignmentSet();
        AlignmentSet(const std::vector<MetaDataTable>& mdts,
                     int fc, int s, int k0, int k1);

            int mc, fc, s, sh, k0, k1, accPix;

            // micrograph < particle < frame <pixels> > >
            std::vector<std::vector<std::vector< Image<T> >>> CCs;
            // micrograph < particle < frame <pixels> > >
            std::vector<std::vector<std::vector< std::vector<gravis::t2Vector<T>> >>> obs;
            // micrograph < particle <pixels> >
            std::vector<std::vector< std::vector<gravis::t2Vector<T>> >> pred;
            // frame <pixels>
            std::vector< std::vector<double> > damage;

            std::vector<std::vector<gravis::d2Vector>> positions;
            std::vector<std::vector<std::vector<gravis::d2Vector>>> initialTracks;
            std::vector<std::vector<gravis::d2Vector>> globComp;

            std::vector<gravis::t2Vector<int>> accCoords;


		template<class T2>
        void copyCC(int m, int p, int f, const Image<T2>& src);
		
        void accelerate(const Image<Complex>& img, std::vector<gravis::t2Vector<T>>& dest);
        void accelerate(const Image<RFLOAT>& img, std::vector<double>& dest);

        gravis::d3Vector updateTsc(
            const std::vector<std::vector<gravis::d2Vector>>& tracks,
            int mg, int threads);
};

template<class T>
AlignmentSet<T>::AlignmentSet()
:   mc(0), fc(0), s(0), sh(0), k0(0), k1(0)
{
}

template<class T>
AlignmentSet<T>::AlignmentSet(
        const std::vector<MetaDataTable> &mdts,
        int fc, int s, int k0, int k1)
:   mc(mdts.size()),
    fc(fc),
    s(s),
    sh(s/2+1),
    k0(k0),
    k1(k1)
{
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
            accCoords.push_back(gravis::t2Vector<int>(x,y));
            num++;
        }
    }

    accPix = num;

    CCs.resize(mc);
    obs.resize(mc);
    pred.resize(mc);

    positions.resize(mc);
    initialTracks.resize(mc);
    globComp.resize(mc);

    for (int m = 0; m < mc; m++)
    {
        const int pc = mdts[m].numberOfObjects();

        positions[m].resize(pc);
        globComp[m].resize(fc);

        initialTracks[m].resize(pc);
        CCs[m].resize(pc);
        obs[m].resize(pc);
        pred[m].resize(pc);

        for (int p = 0; p < pc; p++)
        {
            initialTracks[m][p].resize(fc);
            pred[m][p].resize(accPix);

            CCs[m][p].resize(fc);
            obs[m][p].resize(fc);

            for (int f = 0; f < fc; f++)
            {
                obs[m][p][f].resize(accPix);
            }
        }
    }

    damage.resize(fc);

    for (int f = 0; f < fc; f++)
    {
        damage[f].resize(accPix);
    }
}

template<class T>
template<class T2>
void AlignmentSet<T>::copyCC(int m, int p, int f, const Image<T2> &src)
{
    if (m < 0 || m >= mc ||
        p < 0 || p >= CCs[m].size() ||
        f < 0 || f >= fc)
    {
        REPORT_ERROR_STR("AlignmentSet::copyCC: bad CC-index: "
            << m << ", " << p << ", " << f << " for "
            << mc << ", " << ((m >= 0 && m < mc)? CCs[m].size() : 0) << ", " << fc << ".");
    }

	CCs[m][p][f] = Image<T>(src.data.xdim, src.data.ydim);
	
    for (int y = 0; y < src.data.ydim; y++)
    for (int x = 0; x < src.data.xdim; x++)
    {
        CCs[m][p][f](y,x) = (T)src(y,x);
    }
}

template<class T>
void AlignmentSet<T>::accelerate(
        const Image<Complex> &img,
        std::vector<gravis::t2Vector<T>>& dest)
{
    for (int i = 0; i < accPix; i++)
    {
        gravis::t2Vector<int> c = accCoords[i];

        Complex z = img(c.y, c.x);

        dest[i] = gravis::t2Vector<T>(z.real, z.imag);
    }
}

template<class T>
void AlignmentSet<T>::accelerate(const Image<RFLOAT> &img, std::vector<double>& dest)
{
    for (int i = 0; i < accPix; i++)
    {
        gravis::t2Vector<int> c = accCoords[i];
        dest[i] = img(c.y, c.x);
    }
}

template<class T>
gravis::d3Vector AlignmentSet<T>::updateTsc(
    const std::vector<std::vector<gravis::d2Vector>>& tracks,
    int mg, int threads)
{
    const int pad = 512;
    std::vector<gravis::d3Vector> outT(pad*threads, gravis::d3Vector(0.0, 0.0, 0.0));

    const int pc = tracks.size();

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        int t = omp_get_thread_num();

        const gravis::d2Vector shift = tracks[p][f] / s;

        for (int i = 0; i < accPix; i++)
        {
            gravis::t2Vector<int> acc = accCoords[i];

            double x = acc.x;
            double y = acc.y < sh? acc.y : acc.y - s;

            const double dotp = 2 * PI * (x * shift.x + y * shift.y);

            double a, b;

            SINCOS(dotp, &b, &a);

            const gravis::t2Vector<T> z_obs_t2 = obs[mg][p][f][i];
            const double c = (double) z_obs_t2.x;
            const double d = (double) z_obs_t2.y;

            const double ac = a * c;
            const double bd = b * d;

            const double ab_cd = (a + b) * (c + d);

            const dComplex z_obs(ac - bd, ab_cd - ac - bd);

            const gravis::t2Vector<T> z_pred_t2 = pred[mg][p][i];
            const dComplex z_pred((double) z_pred_t2.x, (double) z_pred_t2.y);

            const double dmg = damage[f][i];

            outT[pad*t][0] += dmg * (z_pred.real * z_obs.real + z_pred.imag * z_obs.imag);
            outT[pad*t][1] += dmg * z_obs.norm();
            outT[pad*t][2] += dmg * z_pred.norm();
        }
    }

    gravis::d3Vector out(0.0, 0.0, 0.0);

    for (int t = 0; t < threads; t++)
    {
        out += outT[pad*t];
    }

    return out;
}



#endif
