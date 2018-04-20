#include <src/jaz/projection_helper.h>
#include <fstream>
#include <cmath>
#include <src/error.h>
#include <string.h>

using namespace gravis;

std::vector<gravis::d4Matrix> ProjectionHelper::loadTiltProjections(
        std::string tiltFile, double centerX, double centerY)
{
    std::ifstream anglesFile(tiltFile.c_str());

    if (!anglesFile.is_open())
    {
        REPORT_ERROR("failed to open " + tiltFile + '\n');
    }

    std::vector<double> angles;
    std::vector<d4Matrix> vol2img;

    const double deg2rad = 3.14159265358979323846/180.0;

    d4Matrix w2i0;
    w2i0(0,3) = -centerX;
    w2i0(1,3) = -centerY;

    while (anglesFile.good())
    {
        double a;
        anglesFile >> a;
        a *= deg2rad;

        angles.push_back(a);

        d4Matrix w2i;
        w2i(0,0) =  cos(a);
        w2i(0,2) =  sin(a);
        w2i(2,0) = -sin(a);
        w2i(2,2) =  cos(a);
        w2i(0,3) =  centerX;
        w2i(1,3) =  centerY;

        vol2img.push_back(w2i*w2i0);
    }

    return vol2img;
}

std::vector<gravis::d4Matrix> ProjectionHelper::loadTiltProjectionsVol(
        std::string tiltFile, double centerX, double centerY,
        double X0, double Y0, double Z0,
        double spacing)
{
    std::ifstream anglesFile(tiltFile.c_str());

    if (!anglesFile.is_open())
    {
        REPORT_ERROR("failed to open " + tiltFile + '\n');
    }

    // vol2world * (0,0,0,1)       = (xw0, yw0, zw0, 1)
    // vol2world * (xwc,ywc,zwc,1) = (xw1, yw1, zw1, 1)

    d4Matrix vol2world;

    vol2world(0,0) = spacing;
    vol2world(1,1) = spacing;
    vol2world(2,2) = spacing;
    vol2world(0,3) = X0;
    vol2world(1,3) = Y0;
    vol2world(2,3) = Z0;

    std::vector<double> angles;
    std::vector<d4Matrix> vol2img;

    const double deg2rad = 3.14159265358979323846/180.0;

    while (anglesFile.good())
    {
        double a;
        anglesFile >> a;
        a *= deg2rad;

        angles.push_back(a);

        d4Matrix w2i;
        w2i(0,0) =  cos(a);
        w2i(0,2) =  sin(a);
        w2i(2,0) = -sin(a);
        w2i(2,2) =  cos(a);
        w2i(0,3) =  centerX;
        w2i(1,3) =  centerY;

        vol2img.push_back(w2i * vol2world);
    }

    return vol2img;
}

std::vector<d4Matrix> ProjectionHelper::loadAffineTransforms(std::string xformFile, double cx, double cy, bool square_result)
{
    std::cout << "img. center: " << cx << ", " << cy << '\n';

    std::ifstream file(xformFile.c_str());

    if (!file.is_open())
    {
        REPORT_ERROR("failed to open " + xformFile + '\n');
    }

    std::vector<d4Matrix> xforms;

    d4Matrix P, Q;
    P.loadIdentity();
    Q.loadIdentity();

    P(0,3) = -cx;
    P(1,3) = -cy;
    Q(0,3) =  cx;

    if (square_result)
    {
        Q(1,3) =  cx;
    }
    else
    {
        Q(1,3) =  cy;
    }

    char text[4096];

    while (file.good())
    {
        file.getline(text, 4096);

        if (strlen(text) < 11) break;

        std::stringstream line(text);

        d4Matrix A;
        A.loadIdentity();

        line >> A(0,0);
        line >> A(0,1);
        line >> A(1,0);
        line >> A(1,1);
        line >> A(0,3);
        line >> A(1,3);

        xforms.push_back(Q*A*P);
    }

    return xforms;
}
