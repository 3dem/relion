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
/*
 Based on rwIMAGIC.h
 Header file for reading and writing Image Science's Imagic files
 Format: 2D image file format for the program Imagic (Image Science)
 Author: Bernard Heymann
 Created: 19990424  Modified: 20011030
*/

#ifndef RWIMAGIC_H_
#define RWIMAGIC_H_

#include "src/metadata_label.h"

#define IMAGICSIZE 1024 // Size of the IMAGIC header for each image

///@defgroup Imagic Imagic File format
///@ingroup ImageFormats

/** Imagic Header
  * @ingroup Imagic
*/
struct IMAGIChead
{             // file header for IMAGIC data
    int imn;          //  0      image location number (1,2,...)
    int ifn;          //  1      # images following
    int ierror;       //  2      error code: error if >0
    int nhfr;         //  3      # header records per image
    int ndate;        //  4      creation day
    int nmonth;       //  5      creation month
    int nyear;        //  6      creation year
    int nhour;        //  7      creation hour
    int nminut;       //  8      creation minute
    int nsec;         //  9      creation second
    int npix2;        // 10      # 4-byte reals in image
    int npixel;       // 11      # image elements
    int ixlp;       // 12      lines per image (Y)
    int iylp;        // 13      pixels per line (X)
    char type[4];      // 14      image type
    int ixold;       // 15      top-left X coordinate
    int iyold;       // 16      top-left Y coordinate
    float avdens;       // 17      average
    float sigma;       // 18      standard deviation
    float varian;       // 19      variance
    float oldavd;      // 20      old average
    float densmax;       // 21      maximum
    float densmin;       // 22      minimum
    //     RFLOAT sum;       // 23+24  sum of densities
    //     RFLOAT squares;    // 25+26  sum of squares
    float dummy[4];   // 23-26  dummy place holder
    char lastpr[8];      // 27+28     last program writing file
    char name[80];       // 29-48     image name
    float extra_1[8];   // 49-56     additional parameters
    float eman_alt;   // 57      EMAN: equiv to psi & PFT omega
    float eman_az;    // 58      EMAN: equiv to theta
    float eman_phi;   // 59      EMAN: equiv to phi
    float extra_2[69];   // 60-128     additional parameters
    float euler_alpha;  // 129   Euler angles: psi
    float euler_beta;  // 130       theta
    float euler_gamma;  // 131       phi
    float proj_weight;  // 132   weight of each projection
    float extra_3[66];   // 133-198     additional parameters
    char history[228];      // 199-255   history
} ;

/************************************************************************
@Function: readIMAGIC
@Description:
 Reading an IMAGIC image format.
@Algorithm:
 A 2D file format for the IMAGIC package.
 The header is stored in a separate file with extension ".hed" and
  a fixed size of 1024 bytes per image.
 The image data is stored in a single block in a file with the
  extension ".img".
 Byte order determination: Year and hour values
        must be less than 256*256.
 Data types:     PACK = byte, INTG = short, REAL = float,
        RECO,COMP = complex float.
 Note that the x and y dimensions are interchanged (actually a display issue).
@Arguments:
 Bimage* p   the image structure.
 int select   image selection in multi-image file (-1 = all images).
@Returns:
 int     error code (<0 means failure).
**************************************************************************/
/** Imagic reader
  * @ingroup Imagic
*/
int  readIMAGIC(long int img_select)
{
#ifdef DEBUG
    printf("DEBUG readIMAGIC: Reading Imagic file\n");
#endif

    IMAGIChead* header = new IMAGIChead;

    if ( fread( header, IMAGICSIZE, 1, fhed ) < 1 )
        REPORT_ERROR((std::string)"readIMAGIC: header file of " + filename + " cannot be read");

    // Determine byte order and swap bytes if from little-endian machine
    char*   b = (char *) header;
    long int extent = IMAGICSIZE - 916;  // exclude char bytes from swapping
    if ( ( abs(header->nyear) > SWAPTRIG ) || ( header->ixlp > SWAPTRIG ) )
    {
        for ( i=0; i<extent; i+=4 )
            if ( i != 56 )          // exclude type string
                swapbytes(b+i, 4);
    }
    long int _xDim,_yDim,_zDim;
    long int _nDim;
    _xDim = (long int) header->iylp;
    _yDim = (long int) header->ixlp;
    _zDim = (long int) 1;
    _nDim = (long int) header->ifn + 1 ;

    std::stringstream Num;
    std::stringstream Num2;
    if ( img_select > (long int)_nDim )
    {
        Num  << img_select;
        Num2 << _nDim;
        REPORT_ERROR((std::string)"readImagic: Image number " + Num.str() +
                     " exceeds stack size " + Num2.str());
    }

    if( img_select > -1)
        _nDim=1;
    data.setDimensions( //setDimensions do not allocate data
        _xDim,
        _yDim,
        _zDim,
        _nDim );
    replaceNsize=_nDim;
    DataType datatype;

    if ( strstr(header->type,"PACK") )
        datatype = UChar;
    else if ( strstr(header->type,"INTG") )
        datatype = SShort;
    else if ( strstr(header->type,"REAL") )
        datatype = Float;
    else if ( strstr(header->type,"RECO") || strstr(header->type,"COMP") )
    {
        REPORT_ERROR("readIMAGIC: only real-space images can be read into RELION");
    }

    // Set min-max values and other statistical values
    if ( header->sigma == 0 && header->varian != 0 )
        header->sigma = sqrt(header->varian);
    if ( header->densmax == 0 && header->densmin == 0 && header->sigma != 0 )
    {
        header->densmin = header->avdens - header->sigma;
        header->densmax = header->avdens + header->sigma;
    }

    MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN,(RFLOAT)header->densmin);
    MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX,(RFLOAT)header->densmax);
    MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG,(RFLOAT)header->avdens);
    MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV,(RFLOAT)header->sigma);
    setSamplingRateInHeader((RFLOAT)1.);
    MDMainHeader.setValue(EMDL_IMAGE_DATATYPE,(int)datatype);

    offset = 0;   // separate header file

    if (dataflag<0)   // Don't read the individual header and the data if not necessary
    {
    	delete header;
    	return 0;
    }

   // Get the header information
    int error_fseek;
    if ( img_select > -1 )
    	error_fseek = fseek( fhed, img_select * IMAGICSIZE, SEEK_SET );
    else
    	error_fseek = fseek( fhed, 0, SEEK_SET );

    if (error_fseek != 0)
    	return -1;

    delete header;

    int pad=0;
    return readData(fimg, img_select, datatype, pad );

}

// The function remained unfinished and broken since 2015. So Takanori removed it in 2021.
void  writeIMAGIC(long int img_select=-1, int mode=WRITE_OVERWRITE, const DataType datatype=Unknown_Type)
{
    REPORT_ERROR("writeIMAGIC() is not implemented.");
}

#endif /* RWIMAGIC_H_ */
