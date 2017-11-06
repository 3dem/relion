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
int readTIFF(TIFF* ftiff, long int img_select, bool isStack=false, const FileName &name="")
{
#undef DEBUG
//#define DEBUG
#ifdef DEBUG
    printf("DEBUG readTIFF: Reading TIFF file. img_select %d\n", img_select);
#endif

    long int _xDim,_yDim,_zDim;
    long int _nDim;

    // libtiff's types
    tsize_t stripSize;
    tstrip_t numberOfStrips;
    uint32 rowsPerStrip, width, length;
    uint16 sampleFormat, bitsPerSample;
    tdata_t buf;
    
    TIFFGetField(ftiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(ftiff, TIFFTAG_IMAGELENGTH, &length);
    _xDim = width;
    _yDim = length;
    _nDim = 1;
    TIFFGetField(ftiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(ftiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetField(ftiff, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);
    stripSize = TIFFStripSize(ftiff);
    numberOfStrips = TIFFNumberOfStrips(ftiff);
    buf = _TIFFmalloc(stripSize);

    // Find the number of frames
    while (TIFFSetDirectory(ftiff, _nDim) != 0) _nDim++;
    //  and back to the start
    TIFFSetDirectory(ftiff, 0);

#ifdef DEBUG
    printf("TIFF width %d, length %d, nDim %d, number of strips %d, strip size %d, sample format %d, bits per sample %d\n", 
           width, length, _nDim, numberOfStrips, stripSize, sampleFormat, bitsPerSample);
#endif

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
        REPORT_ERROR("Unsupport TIFF format\n");
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

    if (img_select != -1) TIFFSetDirectory(ftiff, img_select); // img_select starts from 0
    size_t readsize_n = stripSize * 8 / bitsPerSample;
    size_t haveread_n = 0;
    for (tstrip_t strip = 0; strip < numberOfStrips; strip++) {
        TIFFReadEncodedStrip(ftiff, strip, buf, stripSize);
        castPage2T((char*)buf, MULTIDIM_ARRAY(data) + haveread_n, datatype, readsize_n);
        haveread_n += readsize_n;
    }

    _TIFFfree(buf);

    return 0;
}

#endif
