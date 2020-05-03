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

#include <src/jaz/single_particle/motion_em.h>
#include <src/jaz/single_particle/motion/motion_refinement.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/noise_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/interpolation.h>
#include <omp.h>

using namespace gravis;

MotionEM::MotionEM(
    Projector& projector0,
    Projector& projector1,
    const ObservationModel& obsModel,
    MetaDataTable& viewParams,
    const std::vector<std::vector<Image<Complex>>>& movie,
    const std::vector<gravis::d2Vector>& globalPositions,
    const std::vector<double>& sigma2,
    const std::vector<Image<RFLOAT>>& damageWeights,
    double sig_pos,
    const std::vector<double> &sig_vel_initial,
    const std::vector<double> &sig_div_initial,
    double sig_cutoff,
    int threads)
:
    projector0(projector0),
    projector1(projector1),
    obsModel(obsModel),
    viewParams(viewParams),
    movie(movie),
    globalPositions(globalPositions),
    sigma2(sigma2),
    damageWeights(damageWeights),
    sig_pos(sig_pos),
    sig_vel(movie[0].size() - 1),
    sig_div(movie[0].size() - 1),
    sig_cutoff(sig_cutoff),
    threads(threads),

    fts_full(threads),
    fts_pos(threads),
    fts_vel(threads * sig_vel_initial.size()),
    pc(movie.size()),
    fc(movie[0].size()),
    s_full(movie[0][0].data.ydim),
    sh_full(movie[0][0].data.ydim/2 + 1),
    s_pos(2 * (int)(sig_cutoff * sig_pos)),
    sh_pos((int)(sig_cutoff * sig_pos) + 1),
    s_vel(movie[0].size()-1),
    sh_vel(movie[0].size()-1),
    sig_vel_class(movie[0].size()-1),
    initialized(false)
{
    if (s_pos > s_full)
    {
        s_pos = s_full;
        sh_pos = sh_full;
    }

    for (int f = 0; f < fc-1; f++)
    {
        if (f < sig_vel_initial.size()-1)
        {
            sig_vel_class[f] = f;
            sig_vel[f] = sig_vel_initial[f];
        }
        else
        {
            sig_vel_class[f] = sig_vel_initial.size()-1;
            sig_vel[f] = sig_vel_initial[sig_vel_initial.size()-1];
        }

        if (f < sig_div_initial.size()-1)
        {
            sig_div[f] = sig_div_initial[f];
        }
        else
        {
            sig_div[f] = sig_div_initial[sig_div_initial.size()-1];
        }

        s_vel[f] = 2 * (int)(sig_cutoff * sig_vel[f]);
        sh_vel[f] = (int)(sig_cutoff * sig_vel[f]) + 1;

        if (s_vel[f] > s_pos)
        {
            s_vel[f] = s_pos;
            sh_vel[f] = sh_pos;
        }
    }
}

void MotionEM::computeInitial()
{
    posProb = std::vector<std::vector<Image<RFLOAT>>>(pc);
    velProb = std::vector<std::vector<Image<RFLOAT>>>(pc);
    pred = std::vector<Image<Complex>>(pc);
    initialCC = std::vector<std::vector<Image<RFLOAT>>>(pc);

    for (int p = 0; p < pc; p++)
    {
        initialCC[p] = std::vector<Image<RFLOAT>>(fc, Image<RFLOAT>(s_full,s_full));
        posProb[p] = std::vector<Image<RFLOAT>>(fc);
        velProb[p] = std::vector<Image<RFLOAT>>(fc-1);

        int randSubset;
        viewParams.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
        randSubset -= 1;

        if (randSubset == 0)
        {
            pred[p] = obsModel.predictObservation(projector0, viewParams, p, true, true);
        }
        else
        {
            pred[p] = obsModel.predictObservation(projector1, viewParams, p, true, true);
        }

        MotionRefinement::noiseNormalize(pred[p], sigma2, pred[p]);


        #pragma omp parallel for num_threads(threads)
        for (int f = 0; f < fc; f++)
        {
            std::stringstream stsf;
            stsf << f;

            int threadnum = omp_get_thread_num();

            Image<Complex> ccFs(sh_full,s_full);
            ccFs.data.xinit = 0;
            ccFs.data.yinit = 0;

            for (int y = 0; y < s_full; y++)
            for (int x = 0; x < sh_full; x++)
            {
                ccFs(y,x) = movie[p][f](y,x) * damageWeights[f](y,x) * pred[p](y,x).conj();
            }

            Image<RFLOAT> ccRs(s_full,s_full);

            fts_full[threadnum].inverseFourierTransform(ccFs(), ccRs());

            for (int y = 0; y < s_full; y++)
            for (int x = 0; x < s_full; x++)
            {
                initialCC[p][f](y,x) = s_full * s_full * ccRs(y,x);
            }

            RFLOAT offCC = FilterHelper::maxValue(ccRs);

            ImageOp::linearCombination(ccRs, offCC, 1.0, -1.0, ccRs);

            posProb[p][f] = FilterHelper::expImg(ccRs, s_full * s_full);
            FilterHelper::GaussianEnvelopeCorner2D(posProb[p][f], sig_pos);
            posProb[p][f] = FilterHelper::cropCorner2D(posProb[p][f], s_pos, s_pos);
        }

    } // particles

    initialized = true;
}

void MotionEM::iterate()
{
    updateVelocities();
    consolidateVelocities();
    //smoothVelocities();
    updatePositions(false);
    updatePositions(true);
}

void MotionEM::updateVelocities()
{
    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        int threadnum = omp_get_thread_num();

        std::vector<Image<Complex>> posProbFs(fc);

        for (int f = 0; f < fc; f++)
        {
            fts_pos[threadnum].FourierTransform(posProb[p][f](), posProbFs[f]());
        }

        for (int f = 0; f < fc-1; f++)
        {
            Image<Complex> velProbLargeFs = Image<Complex>(sh_pos, s_pos);

            for (int y = 0; y < s_pos;  y++)
            for (int x = 0; x < sh_pos; x++)
            {
                DIRECT_A2D_ELEM(velProbLargeFs(), y, x) =
                          DIRECT_A2D_ELEM(posProbFs[f+1](), y, x)
                        * DIRECT_A2D_ELEM(posProbFs[f  ](), y, x).conj();
            }

            Image<RFLOAT> velProbLarge(s_pos,s_pos);

            fts_pos[threadnum].inverseFourierTransform(velProbLargeFs(), velProbLarge());

            velProb[p][f] = FilterHelper::cropCorner2D(velProbLarge, s_vel[f], s_vel[f]);
        }
    }
}

void MotionEM::consolidateVelocities(int maxPc)
{
    const int debug_pc = maxPc > 0? maxPc : pc;

    #pragma omp parallel for num_threads(threads)
    for (int f = 0; f < fc-1; f++)
    {
        int threadnum = omp_get_thread_num();

        std::vector<Image<RFLOAT>> velProbNext(pc);
        std::vector<Image<Complex>> velProbFs(pc);

        for (int p = 0; p < pc; p++)
        {
            velProbNext[p] = Image<RFLOAT>(s_vel[f],s_vel[f]);
            fts_vel[threads*sig_vel_class[f] + threadnum].FourierTransform(velProb[p][f](), velProbFs[p]());
        }

        for (int p = 0; p < debug_pc; p++)
        {
            const double eps = 1e-25;

            for (int y = 0; y < s_vel[f]; y++)
            for (int x = 0; x < s_vel[f]; x++)
            {
                if (velProb[p][f](y,x) > eps)
                {
                    velProbNext[p](y,x) = log(velProb[p][f](y,x));
                }
                else
                {
                    velProbNext[p](y,x) = log(eps);
                }
            }

            Image<RFLOAT> velProbB(s_vel[f],s_vel[f]);

            for (int q = 0; q < pc; q++)
            {
                if (q == p) continue;

                const double dist = (globalPositions[p] - globalPositions[q]).length();
                const double sig_real = sig_div[f] * sqrt(dist);

                velProbB = blurVelocity(velProbFs[q], sig_real, f, threadnum);

                for (int y = 0; y < s_vel[f]; y++)
                for (int x = 0; x < s_vel[f]; x++)
                {
                    if (velProbB(y,x) < eps)
                    {
                        velProbNext[p](y,x) += log(eps);
                    }
                    else
                    {
                        velProbNext[p](y,x) += log(velProbB(y,x));
                    }
                }
            }

            double maxVal = -std::numeric_limits<double>::max();

            for (int y = 0; y < s_vel[f]; y++)
            for (int x = 0; x < s_vel[f]; x++)
            {
                const double yy = y < sh_vel[f]? y : y - s_vel[f];
                const double xx = x < sh_vel[f]? x : x - s_vel[f];

                velProbNext[p](y,x) -= 0.5*(xx*xx + yy*yy)/(sig_vel[f]*sig_vel[f]);

                if (velProbNext[p](y,x) > maxVal)
                {
                    maxVal = velProbNext[p](y,x);
                }
            }

            for (int y = 0; y < s_vel[f]; y++)
            for (int x = 0; x < s_vel[f]; x++)
            {
                double v = velProbNext[p](y,x) - maxVal;

                velProbNext[p](y,x) = exp(v);
            }

            velProbNext[p] = NoiseHelper::normalize(velProbNext[p]);
        }

        for (int p = 0; p < debug_pc; p++)
        {
            for (int y = 0; y < s_vel[f]; y++)
            for (int x = 0; x < s_vel[f]; x++)
            {
                velProb[p][f](y,x) = velProbNext[p](y,x);
            }
        }
    }
}


void MotionEM::smoothVelocities()
{
    const double sigma_acc = 1.0;

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        int threadnum = omp_get_thread_num();

        std::vector<Image<RFLOAT>> velProbNext(fc-1);
        std::vector<Image<Complex>> velProbFs(fc-1);

        for (int f = 0; f < fc-1; f++)
        {
            fts_vel[threads*sig_vel_class[f] + threadnum].FourierTransform(velProb[p][f](), velProbFs[f]());
        }

        for (int f = 0; f < fc-1; f++)
        {
            velProbNext[f] = velProb[p][f];

            if (f > 0)
            {
                Image<RFLOAT> velProbOther = adaptSize(
                            blurVelocity(velProbFs[f-1], sigma_acc, f-1, threadnum), s_vel[f]);

                for (int y = 0; y < s_vel[f]; y++)
                for (int x = 0; x < s_vel[f]; x++)
                {
                    velProbNext[f](y,x) *= velProbOther(y,x);
                }
            }

            if (f < fc-2)
            {
                Image<RFLOAT> velProbOther = adaptSize(
                            blurVelocity(velProbFs[f+1], sigma_acc, f+1, threadnum), s_vel[f]);

                for (int y = 0; y < s_vel[f]; y++)
                for (int x = 0; x < s_vel[f]; x++)
                {
                    velProbNext[f](y,x) *= velProbOther(y,x);
                }
            }
        }

        for (int f = 0; f < fc-1; f++)
        {
            for (int y = 0; y < s_vel[f]; y++)
            for (int x = 0; x < s_vel[f]; x++)
            {
                velProb[p][f](y,x) = velProbNext[f](y,x);
            }
        }
    }
}

Image<RFLOAT> MotionEM::blurVelocity(const Image<Complex> &velProbFs, double sigma, int f, int threadnum)
{
    Image<RFLOAT> velProbB(s_vel[f], s_vel[f]);

    const double sig_freq = s_vel[f]/(2.0 * PI * sigma);
    const double sig2_freq = sig_freq * sig_freq;

    Image<Complex> velProbFs_env = velProbFs;

    for (int y = 0; y < s_vel[f]; y++)
    for (int x = 0; x < sh_vel[f]; x++)
    {
        const double yy = y < sh_vel[f]? y : y - s_vel[f];
        const double xx = x;

        velProbFs_env(y,x) *= exp(-0.5*(xx*xx + yy*yy)/sig2_freq);
    }

    fts_vel[threads*sig_vel_class[f] + threadnum].inverseFourierTransform(velProbFs_env(), velProbB());

    return velProbB;
}

Image<RFLOAT> MotionEM::adaptSize(const Image<RFLOAT> &img, int s)
{
    if (img.data.ydim > s)
    {
        return FilterHelper::cropCorner2D(img, s, s);
    }
    else if (img.data.ydim < s)
    {
        return FilterHelper::padCorner2D(img, s, s);
    }
    else
    {
        return img;
    }
}

void MotionEM::updatePositions(bool backward, int maxPc)
{
    const int debug_pc = maxPc > 0? maxPc : pc;

    const int f0 = backward? fc-1 : 0;
    const int f1 = backward? 0    : fc-1;
    const int df = backward? -1   : 1;

    const double eps = 1e-25;

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < debug_pc; p++)
    {
        int threadnum = omp_get_thread_num();

        Image<RFLOAT>  posProbMapped(s_pos, s_pos);
        Image<Complex> posProbFs(sh_pos, s_pos), velProbLargeFs(sh_pos, s_pos);

        for (int f = f0; f != f1; f += df)
        {
            const int ff = f+df;
            const int fv = backward? ff : f;

            Image<RFLOAT> velProbLarge = FilterHelper::padCorner2D(velProb[p][fv], s_pos, s_pos);

            fts_pos[threadnum].FourierTransform(velProbLarge(), velProbLargeFs());
            fts_pos[threadnum].FourierTransform(posProb[p][f](), posProbFs());

            for (int y = 0; y < s_pos; y++)
            for (int x = 0; x < sh_pos; x++)
            {
                if (backward)
                {
                    posProbFs(y,x) = s_pos * s_pos * posProbFs(y,x) * velProbLargeFs(y,x).conj();
                }
                else
                {
                    posProbFs(y,x) = s_pos * s_pos * posProbFs(y,x) * velProbLargeFs(y,x);
                }
            }

            fts_pos[threadnum].inverseFourierTransform(posProbFs(), posProbMapped());

            double sum = 0.0;

            for (int y = 0; y < s_pos; y++)
            for (int x = 0; x < s_pos; x++)
            {
                if (posProbMapped(y,x) < eps)
                {
                    posProb[p][ff](y,x) *= eps;
                }
                else
                {
                    posProb[p][ff](y,x) *= posProbMapped(y,x);
                }

                sum += posProb[p][ff](y,x);
            }

            if (sum > 0.0)
            {
                for (int y = 0; y < s_pos; y++)
                for (int x = 0; x < s_pos; x++)
                {
                    posProb[p][ff](y,x) /= sum;
                }
            }
        }
    }
}

std::vector<d2Vector> MotionEM::getTrack(int particle)
{
    if (!initialized)
    {
        return std::vector<d2Vector>(fc, d2Vector(0.0, 0.0));
    }

    std::vector<d2Vector> out(fc);

    for (int f = 0; f < fc; f++)
    {
        double maxProb = 0.0;
        int bestX = 0, bestY = 0;

        for (int y = 0; y < s_pos; y++)
        for (int x = 0; x < s_pos; x++)
        {
            double p = posProb[particle][f](y,x);

            if (p > maxProb)
            {
                maxProb = p;
                bestX = x;
                bestY = y;
            }
        }

        d2Vector opt(bestX, bestY);

        if (opt.x >= sh_pos)
        {
            opt.x -= s_pos;
        }

        if (opt.y >= sh_pos)
        {
            opt.y -= s_pos;
        }

        out[f] = opt;
    }

    return out;
}

std::vector<d2Vector> MotionEM::getGlobalTrack()
{
    std::vector<d2Vector> out(fc);
    const double eps = 1e-30;

    e_sum = std::vector<Image<RFLOAT>>(fc);

    for (int f = 0; f < fc; f++)
    {
        e_sum[f] = Image<RFLOAT>(s_full, s_full);
        e_sum[f].data.initZeros();

        for (int p = 0; p < pc; p++)
        {
            for (int y = 0; y < s_full; y++)
            for (int x = 0; x < s_full; x++)
            {
                e_sum[f](y,x) += initialCC[p][f](y,x);
            }
        }

        d2Vector pos = Interpolation::quadraticMaxWrapXY(e_sum[f], eps);

        if (pos.x >= sh_full) pos.x -= s_full;
        if (pos.y >= sh_full) pos.y -= s_full;

        out[f] = pos;
    }

    return out;
}
