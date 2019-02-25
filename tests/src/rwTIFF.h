/***************************************************************************
 *
 * Author: "Takanori Nakane"
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

#ifndef RWTIFF_H
#define RWTIFF_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define SAFESET(var, value) ( var = (std::isfinite(value)) ? value : var )

// I/O prototypes
/** TIFF Reader
  * @ingroup TIFF
*/
int readTIFF(TIFF* ftiff, long int img_select, bool readdata=false, bool isStack=false, const FileName &name="")
{
//#define DEBUG_TIFF
#ifdef DEBUG_TIFF
    printf("DEBUG readTIFF: Reading TIFF file. img_select %d\n", img_select);
#endif

    long int _xDim,_yDim,_zDim;
    long int _nDim;

    // libtiff's types
    uint32 width, length;
    uint16 sampleFormat, bitsPerSample;
    
    if (TIFFGetField(ftiff, TIFFTAG_IMAGEWIDTH, &width) != 1 ||
        TIFFGetField(ftiff, TIFFTAG_IMAGELENGTH, &length) != 1)
    {
        REPORT_ERROR("The input TIFF file does not have the width or height field.");
    }
    _xDim = width;
    _yDim = length;
    _zDim = 1;
    _nDim = 1;
    TIFFGetFieldDefaulted(ftiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetFieldDefaulted(ftiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);

    // Find the number of frames
    while (TIFFSetDirectory(ftiff, _nDim) != 0) _nDim++;
    //  and back to the start
    TIFFSetDirectory(ftiff, 0);

#ifdef DEBUG_TIFF
    printf("TIFF width %d, length %d, nDim %d, sample format %d, bits per sample %d\n", 
           width, length, _nDim, sampleFormat, bitsPerSample);
#endif

    // TODO: TIFF is always a stack, isn't it?
    if(isStack)
    {
        _zDim = 1;
        replaceNsize=_nDim;
        std::stringstream Num;
        std::stringstream Num2;
        if (img_select >= (int)_nDim)
        {
            Num  << (img_select + 1);
            Num2 << _nDim;
            REPORT_ERROR((std::string)"readTIFF: Image number " + Num.str() + " exceeds stack size " + Num2.str() + " of image " + name);
        }
    }
    else
        replaceNsize=0;

    // Map the parameters
    if (isStack && img_select==-1)
        _zDim = 1;
    else if(isStack && img_select!=-1)
        _zDim = _nDim = 1;

    data.setDimensions(_xDim, _yDim, _zDim, _nDim);
    data.coreAllocateReuse();
    
    DataType datatype;
    
    if (bitsPerSample == 8 && sampleFormat == 1) {
        datatype = UChar;
    } else if (bitsPerSample == 16 && sampleFormat == 1) {
        datatype = UShort;
    } else if (bitsPerSample == 16 && sampleFormat == 2) {
        datatype = Short;
    } else if (bitsPerSample == 32 && sampleFormat == 3) {
        datatype = Float;
    } else {
        std::cerr << "Unsupported TIFF format: sample format = " << sampleFormat << ", bits per sample = " << bitsPerSample << std::endl;
        REPORT_ERROR("Unsupported TIFF format.\n");
    }
    
    MDMainHeader.setValue(EMDL_IMAGE_DATATYPE,(int)datatype);

    /*
    if ( header->mx && header->a!=0)//ux
        MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X,(RFLOAT)header->a/header->mx);
    if ( header->my && header->b!=0)//yx
        MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y,(RFLOAT)header->b/header->my);
    if ( header->mz && header->c!=0)//zx
        MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z,(RFLOAT)header->c/header->mz);
    */

    if (readdata) {
        if (img_select == -1) img_select = 0; // img_select starts from 0

        size_t haveread_n = 0;
        for (int i = 0; i < _nDim; i++) {
            TIFFSetDirectory(ftiff, img_select);

            // Make sure image property is consistent for all frames
            uint32 cur_width, cur_length;
            uint16 cur_sampleFormat, cur_bitsPerSample;

            if (TIFFGetField(ftiff, TIFFTAG_IMAGEWIDTH, &cur_width) != 1 ||
                TIFFGetField(ftiff, TIFFTAG_IMAGELENGTH, &cur_length) != 1)
            {
                REPORT_ERROR("The input TIFF file does not have the width or height field.");
            }
            TIFFGetFieldDefaulted(ftiff, TIFFTAG_BITSPERSAMPLE, &cur_bitsPerSample);
            TIFFGetFieldDefaulted(ftiff, TIFFTAG_SAMPLEFORMAT, &cur_sampleFormat);
            if ((cur_width != width) || (cur_length != length) || (cur_bitsPerSample != bitsPerSample) ||
                (cur_sampleFormat != sampleFormat)) {
                REPORT_ERROR("All frames in a TIFF should have same width, height and pixel format.\n");
            }

            tsize_t stripSize = TIFFStripSize(ftiff);
            tstrip_t numberOfStrips = TIFFNumberOfStrips(ftiff);
            tdata_t buf = _TIFFmalloc(stripSize);
#ifdef DEBUG_TIFF
            size_t readsize_n = stripSize * 8 / bitsPerSample;
            std::cout << "TIFF stripSize=" << stripSize << " numberOfStrips=" << numberOfStrips << " readsize_n=" << readsize_n << std::endl;
#endif
            for (tstrip_t strip = 0; strip < numberOfStrips; strip++) {
                tsize_t actually_read = TIFFReadEncodedStrip(ftiff, strip, buf, stripSize);
                tsize_t actually_read_n = actually_read * 8 / bitsPerSample;
#ifdef DEBUG_TIFF
                std::cout << "Reading strip: " << strip << "actually read byte:" << actually_read << std::endl;
#endif
                castPage2T((char*)buf, MULTIDIM_ARRAY(data) + haveread_n, datatype, actually_read_n);
                haveread_n += actually_read_n;
	    }

            _TIFFfree(buf);
            img_select++;
        }

        /* Flip the Y axis.
 
           In an MRC file, the origin is bottom-left, +X to the right, +Y to the top.
           (c.f. Fig. 2 of Heymann et al, JSB 2005 https://doi.org/10.1016/j.jsb.2005.06.001
                 IMOD's interpretation http://bio3d.colorado.edu/imod/doc/mrc_format.txt)
           3dmod (from IMOD) and e2display.py (from EMAN2) display like this.

           relion_display has the origin at top-left, +X to the right, +Y to the bottom.
           GIMP and ImageJ display in this way as well.
           A TIFF file, with TIFFTAG_ORIENTATION = 1 (default), shares this convention.

           So, the origin and the direction of the Y axis are the opposite between MRC and TIFF.
           IMOD, EMAN2, SerialEM and MotionCor2 flip the Y axis whenever they read or write a TIFF file.
           We follow this.
        */

        T tmp;
        const int ylim = _yDim / 2, z = 0;
        for (int n = 0; n < _nDim; n++) {
            for (int y1 = 0; y1 < ylim; y1++) {
                const int y2 = _yDim - 1 - y1;
                for (int x = 0; x < _xDim; x++) { // TODO: memcpy or pointer arithmetic is probably faster
                    tmp = DIRECT_NZYX_ELEM(data, n, z, y1, x);
                    DIRECT_NZYX_ELEM(data, n, z, y1, x) = DIRECT_NZYX_ELEM(data, n, z, y2, x);
                    DIRECT_NZYX_ELEM(data, n, z, y2, x) = tmp;
                }
            }
        } 
    }

    return 0;
}

#endif
