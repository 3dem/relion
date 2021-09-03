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

#ifndef MOTION_ESTIMATOR_H
#define MOTION_ESTIMATOR_H

#include <src/image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class IOParser;
class ParFourierTransformer;
class ReferenceMap;
class ObservationModel;
class MicrographHandler;

class MotionEstimator
{
    public:

        MotionEstimator();


        void read(IOParser& parser, int argc, char *argv[]);

        void init(int verb, int fc, int nr_omp_threads,
                  bool debug, std::string outPath,
                  ReferenceMap* reference,
                  ObservationModel* obsModel,
                  MicrographHandler* micrographHandler);

        void process(const std::vector<MetaDataTable> &mdts, long g_start, long g_end, bool do_update_FCC);


        // load micrograph from mdt and compute all data required for the optimization;
        // positions, initialTracks and globComp need to have the right sizes already (pc, pc*fc, fc)
        // (also used by MotionParamEstimator)
        void prepMicrograph(
            // in:
            const MetaDataTable& mdt, std::vector<ParFourierTransformer>& fts,
			const std::vector<Image<RFLOAT>>& dmgWeight,
			int ogmg,
            // out:
            std::vector<std::vector<Image<Complex>>>& movie,
            std::vector<std::vector<Image<RFLOAT>>>& movieCC,
            std::vector<gravis::d2Vector>& positions,
            std::vector<std::vector<gravis::d2Vector>>& initialTracks,
            std::vector<gravis::d2Vector>& globComp);

        // perform the actual optimization (also used by MotionParamEstimator)
        std::vector<std::vector<gravis::d2Vector>> optimize(
            const std::vector<std::vector<Image<double>>>& movieCC,
            const std::vector<std::vector<gravis::d2Vector>>& inTracks,
            double sig_vel_px, double sig_acc_px, double sig_div_px,
            const std::vector<gravis::d2Vector>& positions,
            const std::vector<gravis::d2Vector>& globComp) const;

        // syntactic sugar for float-valued CCs
        std::vector<std::vector<gravis::d2Vector>> optimize(
            const std::vector<std::vector<Image<float>>>& movieCC,
            const std::vector<std::vector<gravis::d2Vector>>& inTracks,
            double sig_vel_px, double sig_acc_px, double sig_div_px,
            const std::vector<gravis::d2Vector>& positions,
            const std::vector<gravis::d2Vector>& globComp) const;

	std::vector<Image<RFLOAT>> computeDamageWeights(int opticsGroup);
		
        bool isReady();

        double getDosePerFrame();
        void proposeDosePerFrame(double dpf, std::string metaFn, int verb);

		double getCCPad();

		static std::vector<bool> findUnfinishedJobs(
                const std::vector<MetaDataTable>& mdts, std::string path);

        // translates the given parameters (in A or A/dose) into pixels
        // done in one place to ensure consistency
        double normalizeSigVel(double sig_vel, double angpix);
        double normalizeSigDiv(double sig_div, double angpix);
        double normalizeSigAcc(double sig_acc, double angpix);
		
		int getVerbosity();
		void setVerbosity(int v);


    protected:

            bool paramsRead, ready;

            // read from cmd line
            int maxEDs, maxIters, globOffMax, group;

            bool unregGlob, globOff, cutoffOut,
                diag, expKer, global_init, debugOpt,
				params_scaled_by_dose;

            double dmga, dmgb, dmgc, dosePerFrame,
                sig_vel, sig_div, sig_acc, optEps,
				cc_pad;
			
			std::string paramsFn;

			// @TODO: allow for varying fc (frame count)
            // set at init
            int fc, verb, nr_omp_threads, s_ref, sh_ref;
			std::vector<int> s, sh;
			
			double angpix_ref;
            std::vector<double> angpix;			
            bool debug, no_whitening, all_groups;
			
			std::vector<std::vector<Image<RFLOAT>>> damageWeights;

            std::string outPath;

            ReferenceMap* reference;
            ObservationModel* obsModel;
            MicrographHandler* micrographHandler;


        void updateFCC(
            const std::vector<std::vector<Image<Complex>>>& movie,
            const std::vector<std::vector<gravis::d2Vector>>& tracks,
            const MetaDataTable& mdt,
            std::vector<Image<RFLOAT>>& tables,
            std::vector<Image<RFLOAT>>& weights0,
            std::vector<Image<RFLOAT>>& weights1);
		
		void writeFCC(
			const std::vector<Image<RFLOAT>>& fccData,
			const std::vector<Image<RFLOAT>>& fccWeight0,
			const std::vector<Image<RFLOAT>>& fccWeight1,
			std::string fn_root);

		void writeTracks(
            const std::vector<std::vector<gravis::d2Vector>>& tracks,
			double angpix_mg,
            const std::vector<gravis::d2Vector>& positions,
            std::string fn_root, double visScale);


        static bool isJobFinished(std::string filenameRoot);

};

#endif
