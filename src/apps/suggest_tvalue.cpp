/***************************************************************************
*
* Author: "Sjors H.W. Scheres"
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

#include <src/args.h>
#include "src/time.h"
#include <src/image.h>
#include <src/transformations.h>
#include <src/fftw.h>
#include <src/funcs.h>

class SuggestTvalue {
public:
    // I/O Parser
    IOParser parser;

    FileName fn_map, fn_mask;
    RFLOAT standard_t;
    int n_try;
public:
	SuggestTvalue() { }

	// Read command line arguments
	void read(int argc, char **argv);

	// Execute
	void run();

    RFLOAT getSumAmpl2(MultidimArray<RFLOAT> map);

};

void SuggestTvalue::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

    int general_section = parser.addSection("Options");
    fn_map = parser.getOption("--map", "Consensus map");
    fn_mask = parser.getOption("--mask", "Mask used for focussed classification/refinement");
    standard_t	= textToFloat(parser.getOption("--T", "Standard T-value", "4"));
    n_try = textToInteger(parser.getOption("--try", "Number of times to position mask randomly to find a high-power density area", "10"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}
RFLOAT SuggestTvalue::getSumAmpl2(MultidimArray<RFLOAT> map)
{

    FourierTransformer transformer;
    MultidimArray<Complex > FT;
	transformer.FourierTransform(map, FT, false);

    RFLOAT result = 0.;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
            result += norm(DIRECT_A3D_ELEM(FT, k, i, j));
	}

    return result*NZYXSIZE(map);

}
void SuggestTvalue::run()
{

    randomize_random_generator();

    Image<RFLOAT> Imap, Imask;
    Imap.read(fn_map);
    Imask.read(fn_mask);
    Imap().setXmippOrigin();
	Imask().setXmippOrigin();

	if (!Imap().sameShape(Imask()))
	{
		std::cerr << " Size of map: "; Imap().printShape(std::cerr); std::cerr << std::endl;
		std::cerr << " Size of mask: "; Imask().printShape(std::cerr); std::cerr << std::endl;
		REPORT_ERROR("SuggestTvalue ERROR: The map and mask are not of the same size!");
	}

    // Check values are between 0 and 1
    RFLOAT avg, stddev, minval, maxval;
    Imask().computeStats(avg, stddev, minval, maxval);

    if (minval < -1e-6 || maxval - 1. > 1.e-6)
    {
        std::cerr << " minval= " << minval << " maxval= " << maxval << std::endl;
        REPORT_ERROR("SuggestTvalue ERROR: mask values not in range [0,1]!");
    }
    std::cout << "  This program provides suggestions for regularization parameter T in focussed classification or refinement " << std::endl;
    std::cout << "  The mask occupies " << avg*100. << "% of the box" << std::endl;
    std::cout << "  Standard T-value = " << standard_t << std::endl;

	RFLOAT sum_ori = getSumAmpl2(Imap());
	RFLOAT sum_msk = getSumAmpl2(Imap() * Imask());

    std::cout << "  Randomly shifting mask around the map to find a region with higher power. " << std::endl;
    init_progress_bar(n_try+1);
    // Place mask in the center of the map
    selfTranslateCenterOfMassToCenter(Imask());
    RFLOAT sum_max = getSumAmpl2(Imap() * Imask());

    // Translate mask n_try times randomly within 1/3 of the box to see if there are stronger powers than at the center
    MultidimArray<RFLOAT> Mshift;
    Matrix1D<RFLOAT> shift(3);
    RFLOAT thirdbox = XSIZE(Imap())/3.;
    for (int i = 0; i < n_try; i++)
    {
        XX(shift) = ROUND(rnd_unif(-thirdbox, thirdbox));
        YY(shift) = ROUND(rnd_unif(-thirdbox, thirdbox));
        ZZ(shift) = ROUND(rnd_unif(-thirdbox, thirdbox));
        translate(Imask(), Mshift, shift);
        RFLOAT sum_shift = getSumAmpl2(Imap() * Mshift);
        if (sum_shift > sum_max) sum_max = sum_shift;
        //std::cout << " random shift = " << shift << " sum_shift " << sum_shift << std::endl;

        progress_bar(i+1);
    }
    progress_bar(n_try+1);

    std::cout << "  " << std::endl;
    std::cout << "  Power of input map                        = " << sum_ori << std::endl;
    std::cout << "  Power of masked map                       = " << sum_msk << std::endl;
    std::cout << "  Max power of randomly shifted masked map  = " << sum_max << std::endl;

    RFLOAT t_min = XMIPP_MIN(standard_t * sum_ori/sum_max, standard_t * sum_ori/sum_msk);
    RFLOAT t_max = XMIPP_MAX(standard_t * sum_ori/sum_max, standard_t * sum_ori/sum_msk);

    std::cout << "  " << std::endl;
    std::cout << "  Suggested T-value based on power in masked region: " << standard_t * sum_ori/sum_msk << std::endl;
    std::cout << "  Suggested T-value based on power in other regions: " << standard_t * sum_ori/sum_max << std::endl;
    std::cout << "  " << std::endl;
    std::cout << "  You may want to run focussed classifications or refinements with T-values around these values... " << std::endl;
    std::cout << "  But, always keep an eye on the noise in the resulting maps: " << std::endl;
    std::cout << "  Too much high-res noise means your T-value is too high." << std::endl;
    std::cout << "  Too low-resolution maps mean your T-value is too low." << std::endl;


}

int main(int argc, char *argv[])
{

    SuggestTvalue app;

    try {
        app.read(argc, argv);
        app.run();

    }
    catch (RelionError XE) {
        std::cerr << XE;
        return RELION_EXIT_FAILURE;
    }
    return RELION_EXIT_SUCCESS;

}

