#ifndef REFERENCE_MAP_H
#define REFERENCE_MAP_H

#include <src/image.h>
#include <src/projector.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class ObservationModel;
class LegacyObservationModel;

class ReferenceMap
{
    public:

        typedef enum {Own, Opposite} HalfSet;

        ReferenceMap();
	
			// input parameters:
            std::string reconFn0, reconFn1, maskFn, fscFn;
            double paddingFactor;

            // data:
            Image<RFLOAT> freqWeight, mask;
            std::vector<double> freqWeight1D;
			Projector projectors[2];
            int k_out, s, sh;
			bool hasMask;

        void read(IOParser& parser, int argc, char *argv[]);
        void load(int verb, bool debug);
		
		Image<RFLOAT> getHollowWeight(double kmin_px);
		
		std::vector<Image<Complex>> predictAll(
				const MetaDataTable& mdt,
                ObservationModel& obs,
                HalfSet hs, int threads,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);

        Image<Complex> predict(
                const MetaDataTable& mdt, int p,
                ObservationModel& obs,
                HalfSet hs,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);
		
		std::vector<Volume<gravis::t2Vector<Complex> > > predictAllComplexGradients(
				const MetaDataTable& mdt,
                ObservationModel& obs,
                HalfSet hs, int threads,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);
		
		Volume<gravis::t2Vector<Complex>> predictComplexGradient(
                const MetaDataTable& mdt, int p,
                ObservationModel& obs,
                HalfSet hs,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);
		
		std::vector<Image<Complex>> predictAll(
				const MetaDataTable& mdt,
                const LegacyObservationModel& obs,
                HalfSet hs, int threads,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);

        Image<Complex> predict(
                const MetaDataTable& mdt, int p,
                const LegacyObservationModel& obs,
                HalfSet hs,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);
};

#endif
