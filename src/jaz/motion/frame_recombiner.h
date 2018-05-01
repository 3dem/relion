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
                  int verb, int s, int fc, 
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
            int k0, k1;
            double k0a, k1a, b_scale;
            std::string bfacFn;

            // set at init:
            int s, sh, fc;
            int verb, nr_omp_threads;
            std::string outPath;
            bool debug;
            double angpix, maxFreq;

            ObservationModel* obsModel;
            MicrographHandler* micrographHandler;

            // computed by weightsFromFCC or weightsFromBfacs:
            std::vector<Image<RFLOAT>> freqWeights;


        std::vector<Image<RFLOAT>> weightsFromFCC(const std::vector<MetaDataTable>& allMdts);
        std::vector<Image<RFLOAT>> weightsFromBfacs();

        static bool isJobFinished(std::string filenameRoot);
};

#endif
