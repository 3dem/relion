#ifndef ALIGNMENT_SET_H
#define ALIGNMENT_SET_H

#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t3Vector.h>
#include <vector>

class AlignmentSet
{
    public:

        AlignmentSet();
        AlignmentSet(const std::vector<MetaDataTable>& mdts, int fc, int s, int k0, int k1);

            int mc, fc, s, sh, k0, k1, accPix;

            // micrograph < particle < frame <pixels> > >
            std::vector<std::vector<std::vector< Image<double> >>> CCs;            
            // micrograph < particle < frame <pixels> > >
            std::vector<std::vector<std::vector< std::vector<gravis::d2Vector> >>> obs;
            // micrograph < particle <pixels> >
            std::vector<std::vector< std::vector<gravis::d2Vector> >> pred;
            // frame <pixels>
            std::vector< std::vector<double> > damage;

            std::vector<std::vector<gravis::d2Vector>> positions;
            std::vector<std::vector<std::vector<gravis::d2Vector>>> initialTracks;
            std::vector<std::vector<gravis::d2Vector>> globComp;

            std::vector<gravis::t2Vector<int>> accCoords;

        std::vector<gravis::d2Vector> accelerate(const Image<Complex>& img);
        std::vector<double> accelerate(const Image<RFLOAT>& img);

        gravis::d3Vector updateTsc(
            const std::vector<std::vector<gravis::d2Vector>>& tracks,
            int mg, int threads);
};


#endif
