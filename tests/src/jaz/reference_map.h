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

#ifndef REFERENCE_MAP_H
#define REFERENCE_MAP_H

#include <src/image.h>
#include <src/projector.h>
#include <vector>

class ObservationModel;

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
                const ObservationModel& obs,
                HalfSet hs, int threads,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);

        Image<Complex> predict(
                const MetaDataTable& mdt, int p,
                const ObservationModel& obs,
                HalfSet hs,
                bool applyCtf = true,
                bool applyTilt = true, 
				bool applyShift = true);

};

#endif
