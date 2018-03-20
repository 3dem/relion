#include <src/jaz/ctf_helper.h>
#include <src/jaz/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>


std::vector<CTF> CtfHelper :: loadCtffind4(std::string path, int imageCount, double voltage, double Cs, double Q0, double Bfac, double scale)
{
    /*
     example:
    # Output from CTFFind version 4.1.5, run on 2017-03-30 15:12:45
    # Input file: /beegfs/zivanov/tomograms/ts_05/frames/05_f32.mrc ; Number of micrographs: 1
    # Pixel size: 1.000 Angstroms ; acceleration voltage: 300.0 keV ; spherical aberration: 2.70 mm ; amplitude contrast: 0.07
    # Box size: 512 pixels ; min. res.: 30.0 Angstroms ; max. res.: 5.0 Angstroms ; min. def.: 5000.0 um; max. def. 50000.0 um
    # Columns: #1 - micrograph number; #2 - defocus 1 [Angstroms]; #3 - defocus 2; #4 - azimuth of astigmatism; #5 - additional phase shift [radians]; #6 - cross correlation; #7 - spacing (in Angstroms) up to which CTF rings were fit successfully
    1.000000 10295.926758 10012.275391 -38.856349 0.000000 0.030650 5.279412
    */

    size_t ast = path.find_first_of('*');
    if (ast == std::string::npos)
    {
        REPORT_ERROR("CtfHelper::loadCtffind4: asterisk required in path.\n");
    }

    std::string fnBase = path.substr(0, ast);
    std::string fnEnd = path.substr(ast+1);

    std::vector<CTF> ctfs(imageCount);

    char text[4096];

    for (int i = 0; i < imageCount; i++)
    {
        std::stringstream sts;
        sts << i;
        std::string fnm;
        sts >> fnm;

        std::string fn = fnBase + fnm + fnEnd;
        std::ifstream file(fn.c_str());

        if (!file.is_open())
        {
            REPORT_ERROR("failed to open " + fn + '\n');
        }

        while (file.getline(text, 4096))
        {
            if (text[0] == '#') continue;

            std::stringstream line(text);

            /*
            #1 - micrograph number;
            #2 - defocus 1 [Angstroms];
            #3 - defocus 2;
            #4 - azimuth of astigmatism;
            #5 - additional phase shift [radians];
            #6 - cross correlation;
            #7 - spacing (in Angstroms) up to which CTF rings were fit successfully
            */

            double imgNumber, defocus1, defocus2, azimuth, phaseShift, crossCorr, bestBefore;

            line >> imgNumber;
            line >> defocus1;
            line >> defocus2;
            line >> azimuth;
            line >> phaseShift;
            line >> crossCorr;
            line >> bestBefore;

            ctfs[i].setValues(defocus1, defocus2, azimuth, voltage, Cs, Q0, Bfac, scale, phaseShift);

                    /*RFLOAT _defU, RFLOAT _defV, RFLOAT _defAng,
                    RFLOAT _voltage, RFLOAT _Cs, RFLOAT _Q0,
                    RFLOAT _Bfac, RFLOAT _scale = 1., RFLOAT _phase_shift = 0.*/
        }
    }

    return ctfs;
}
