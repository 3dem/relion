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

#ifndef FRAME_RECOMBINER_H
#define FRAME_RECOMBINER_H

#include <src/image.h>
#include <vector>
#include <string>

class IOParser;
class ObservationModel;
class MicrographHandler;

class FrameRecombiner
{
    public:

        FrameRecombiner();


        void read(IOParser& parser, int argc, char *argv[]);

        void init(const std::vector<MetaDataTable>& allMdts,
                  int verb, int s_ref, int fc, 
				  double maxFreq, int nr_omp_threads,
                  std::string outPath, bool debug,
                  ObservationModel* obsModel,
                  MicrographHandler* micrographHandler);

        void process(const std::vector<MetaDataTable>& mdts, long g_start, long g_end);


        bool doingRecombination();
		
		// has a max. freq. parameter been supplied?
		bool outerFreqKnown();


        static std::vector<MetaDataTable> findUnfinishedJobs(
                const std::vector<MetaDataTable>& mdts, std::string path);


    protected:

            // read from cmd. line:
            bool doCombineFrames, bfac_diag;
            int k0, k1, box_arg, scale_arg;
            double k0a, k1a;
            std::string bfacFn;

            // set at init:
            int s_ref, sh_ref, fc, s_out, sh_out;
            int verb, nr_omp_threads;
            std::string outPath;
            bool debug;
            double angpix_ref, angpix_out, maxFreq;

            ObservationModel* obsModel;
            MicrographHandler* micrographHandler;

            // computed by weightsFromFCC or weightsFromBfacs:
            std::vector<Image<RFLOAT>> freqWeights;


        std::vector<Image<RFLOAT>> weightsFromFCC(const std::vector<MetaDataTable>& allMdts);
        std::vector<Image<RFLOAT>> weightsFromBfacs();

        static bool isJobFinished(std::string filenameRoot);
};

#endif
