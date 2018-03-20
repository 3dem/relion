#ifndef COMPLEX_IO_H
#define COMPLEX_IO_H

#include <src/image.h>
#include <src/complex.h>
#include <string>

class ComplexIO
{
    public:

        static void write(const MultidimArray<Complex>& img, std::string fnBase, std::string fnSuffix);
        static void read(Image<Complex>& img, std::string fnBase, std::string fnSuffix);

};

#endif
