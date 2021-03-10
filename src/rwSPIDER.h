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
 Based on rwSPIDER.h
 Header file for reading and writing SPIDER files
 Format: 3D image file format for the SPIDER package
 Author: Bernard Heymann
 Created: 19990410  Modified: 20010928
*/

#ifndef RWSPIDER_H
#define RWSPIDER_H

#define SPIDERSIZE 1024 // Minimum size of the SPIDER header (variable)
///@defgroup Spider Spider File format
///@ingroup ImageFormats

/** Spider Header
  * @ingroup Spider
*/
struct SPIDERhead
{                    // file header for SPIDER data
    float nslice;    //  0      slices in volume (image = 1)
    float nrow;      //  1      rows per slice
    float irec;      //  2      # records in file (unused)
    float nhistrec;  //  3      (obsolete)
    float iform;     //  4      file type specifier
    float imami;     //  5      max/min flag (=1 if calculated)
    float fmax;      //  6      maximum
    float fmin;      //  7      minimum
    float av;        //  8      average
    float sig;       //  9      standard deviation (=-1 if not calculated)
    float ihist;     // 10      (obsolete)
    float nsam;      // 11      pixels per row
    float labrec;    // 12      # records in header
    float iangle;    // 13      flag: tilt angles filled
    float phi;       // 14      tilt angles
    float theta;     // 15
    float gamma;     // 16      (=psi)
    float xoff;      // 17      translation
    float yoff;      // 18
    float zoff;      // 19
    float scale;     // 20      scaling
    float labbyt;    // 21      # bytes in header
    float lenbyt;    // 22      record length in bytes (row length)
    float istack;    // 23      indicates stack of images
    float inuse;     // 24      indicates this image in stack is used (not used)
    float maxim;     // 25      max image in stack used
    float imgnum;    // 26      number of current image
    float unused[2]; // 27-28     (unused)
    float kangle;    // 29      flag: additional angles set
    float phi1;      // 30      additional angles
    float theta1;    // 31
    float psi1;      // 32
    float phi2;      // 33
    float theta2;    // 34
    float psi2;      // 35

    double fGeo_matrix[3][3]; // x9 = 72 bytes: Geometric info
    float fAngle1; // angle info

    float fr1;
    float fr2; // lift up cosine mask parameters

    /** Fraga 23/05/97  For Radon transforms **/
    float RTflag; // 1=RT, 2=FFT(RT)
    float Astart;
    float Aend;
    float Ainc;
    float Rsigma; // 4*7 = 28 bytes
    float Tstart;
    float Tend;
    float Tinc; // 4*3 = 12, 12+28 = 40B

    char fNada2[576]; // empty 700-76-40=624-40-8= 576 bytes

    char cdat[12];   // 211-213   creation date
    char ctim[9];  // 214-215   creation time
    char ctit[160];  // 216-255   title
} ;

/************************************************************************
@Function: readSPIDER
@Description:
 Reading a SPIDER image file format.
@Algorithm:
 A 3D multi-image format used in electron microscopy.
 Header size:    1024 bytes (not same as data offset!).
 Data offset:    sizeof(float)*x_size*ceil(1024/x_size)
 File format extensions:   .spi
 Byte order determination: File type and third dimension values
        must be less than 256*256.
 Data type:      only float.
 A multi-image file has a global header followed by a header and data
 for each sub-image.
@Arguments:
 Bimage* p   the image structure.
 int img_select  image selection in multi-image file (-1 = all images).
@Returns:
 int     error code (<0 means failure).
**************************************************************************/
/** Spider Reader
  * @ingroup Spider
*/

int readSPIDER(long int img_select)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readSPIDER: Reading Spider file\n");
#endif
#undef DEBUG

    SPIDERhead* header = new SPIDERhead;
    if ( fread( header, SPIDERSIZE, 1, fimg ) < 1 )
        REPORT_ERROR("rwSPIDER: cannot allocate memory for header");

    swap = 0;

    // Determine byte order and swap bytes if from different-endian machine
    char*    b = (char *) header;
    int      i;
    int      extent = SPIDERSIZE - 180;  // exclude char bytes from swapping
    if ( ( fabs(header->nrow) > SWAPTRIG ) || ( fabs(header->iform) > SWAPTRIG ) ||
         ( fabs(header->nslice) < 1 ) )
    {
        swap = 1;
        for ( i=0; i<extent; i+=4 )
            swapbytes(b+i, 4);
    }

    if(header->labbyt != header->labrec*header->lenbyt)
        REPORT_ERROR((std::string)"Invalid Spider file:  " + filename);

    offset = (int) header->labbyt;
    DataType datatype  = Float;

    MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN,(RFLOAT)header->fmin);
    MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX,(RFLOAT)header->fmax);
    MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG,(RFLOAT)header->av);
    MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV,(RFLOAT)header->sig);
    setSamplingRateInHeader((RFLOAT)header->scale);

    bool isStack = ( header->istack > 0 );
    long int _xDim,_yDim,_zDim;
    long int _nDim, _nDimSet;
    _xDim = (long int) header->nsam;
    _yDim = (long int) header->nrow;
    _zDim = (long int) header->nslice;
    _nDim = 1;

    if(isStack)
    {
        _nDim = (long int) header->maxim;
        replaceNsize=_nDim;
    }
    else
        replaceNsize=0;

    /************
     * BELOW HERE DO NOT USE HEADER BUT LOCAL VARIABLES
     */

   // Map the parameters, REad the whole object (-1) or a slide
    // Only handle stacks of images not of volumes
    if(!isStack)
        _nDimSet = 1;
    else
    {
        if(img_select==-1)
            _nDimSet = _nDim;
        else
            _nDimSet = 1;
    }

    data.setDimensions(_xDim, _yDim, _zDim, _nDimSet);

    if (isStack && dataflag<0)   // Don't read the individual header and the data if not necessary
    {
    	delete header;
    	return 0;
    }

    size_t pad         = 0;

    std::stringstream Num;
    std::stringstream Num2;
    //image is in stack? and set right initial and final image
    if ( isStack)
    {
        pad         = offset;
        if ( img_select > _nDim )
        {
            Num  << img_select;
            Num2 << _nDim;
            REPORT_ERROR((std::string)"readSpider: Image number " + Num.str() +
                         " exceeds stack size " + Num2.str());
        }
        offset += offset;
    }

    delete header;

#ifdef DEBUG
    size_t header_size = offset;
    size_t image_size  = header_size + ZYXSIZE(data)*sizeof(float);
    std::cerr<<"DEBUG readSPIDER: header_size = "<<header_size<<" image_size = "<<image_size<<std::endl;
    std::cerr<<"DEBUG readSPIDER: img_select= "<<img_select<<" n= "<<Ndim<<" pad = "<<pad<<std::endl;
#endif
    //offset should point to the begin of the data
    return readData(fimg, img_select, datatype, pad );

}
/************************************************************************
@Function: writeSPIDER
@Description:
 Writing a SPIDER image file format.
@Algorithm:
 A 3D image format used in electron microscopy.
@Arguments:
@Returns:
 int     error code (<0 means failure).
**************************************************************************/
/** Spider Writer
  * @ingroup Spider
*/
int  writeSPIDER(long int select_img=-1, bool isStack=false, int mode=WRITE_OVERWRITE, const DataType dtype=Unknown_Type)
{
#undef DEBUG
//#define DEBUG
#ifdef DEBUG
    printf("DEBUG writeSPIDER: Writing Spider file\n");
    printf("DEBUG writeSPIDER: File %s\n", filename.c_str());
#endif
//#undef DEBUG

    if (dtype != Unknown_Type && dtype != Float)
        REPORT_ERROR("writeSPIDER() can write only in Float.");

    //check if we are going to add or substitute an slice
    //in an existing stack
    //IsStack?
    //else
    long int Xdim = XSIZE(data);
    long int Ydim = YSIZE(data);
    long int Zdim = ZSIZE(data);
    long int Ndim = NSIZE(data);

    float  lenbyt = sizeof(float)*Xdim;  // Record length (in bytes)
    float  labrec = floor(SPIDERSIZE/lenbyt); // # header records
    if ( fmod(SPIDERSIZE,lenbyt) != 0 )
        labrec++;
    float  labbyt = labrec*lenbyt;   // Size of header in bytes
    offset = (int) labbyt;
    SPIDERhead* header = (SPIDERhead *) askMemory((int)labbyt*sizeof(char));

    // Map the parameters
    header->lenbyt = lenbyt;     // Record length (in bytes)
    header->labrec = labrec;     // # header records
    header->labbyt = labbyt;     // Size of header in bytes

    header->irec   = labrec + floor((ZYXSIZE(data)*sizeof(float))/lenbyt + 0.999999); // Total # records
    header->nsam   = Xdim;
    header->nrow   = Ydim;
    header->nslice = Zdim;

#ifdef DEBUG
    printf("DEBUG writeSPIDER: Size: %g %g %g\n", header->nsam, header->nrow, header->nslice);
#endif

    if ( Zdim < 2 )
    	header->iform = 1;     // 2D image
    else
    	header->iform = 3;     // 3D volume
    RFLOAT aux;
    header->imami = 0;//never trust max/min

    if (!MDMainHeader.isEmpty())
    {
#ifdef DEBUG
    	std::cerr<<"Non-empty MDMainHeader"<<std::endl;
#endif
    	if(MDMainHeader.getValue(EMDL_IMAGE_STATS_MIN,   aux))
            header->fmin = (float)aux;
        if(MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX,   aux))
            header->fmax = (float)aux;
        if(MDMainHeader.getValue(EMDL_IMAGE_STATS_AVG,   aux))
            header->av   = (float)aux;
        if(MDMainHeader.getValue(EMDL_IMAGE_STATS_STDDEV,aux))
            header->sig  = (float)aux;
    }
    // For multi-image files
    if (Ndim > 1 || mode == WRITE_APPEND || isStack)
    {
        header->istack = 2;
        header->inuse =  1;
        header->maxim = Ndim;
        if(mode == WRITE_APPEND)
            header->maxim = replaceNsize +1;
    }
    else
    {
        header->istack = 0;
        header->inuse = 0;
        header->maxim = 1;
    }

    //else end
    // Set time and date
    time_t timer;
    time ( &timer );
    tm* t = localtime(&timer);
    while ( t->tm_year > 100 )
        t->tm_year -= 100;
    sprintf(header->ctim, "%02d:%02d:%02d", t->tm_hour, t->tm_min, t->tm_sec);
    sprintf(header->cdat, "%02d-%02d-%02d", t->tm_mday, t->tm_mon, t->tm_year);

    size_t datasize, datasize_n;
    datasize_n = Xdim*Ydim*Zdim;
    datasize = datasize_n * gettypesize(Float);

#ifdef DEBUG

    printf("DEBUG writeSPIDER: Date and time: %s %s\n", header->cdat, header->ctim);
    printf("DEBUG writeSPIDER: Text label: %s\n", header->ctit);
    printf("DEBUG writeSPIDER: Header size: %g\n", header->labbyt);
    printf("DEBUG writeSPIDER: Header records and record length: %g %g\n", header->labrec, header->lenbyt);
    printf("DEBUG writeSPIDER: Data size: %ld\n", datasize);
    printf("DEBUG writeSPIDER: Data offset: %ld\n", offset);
    printf("DEBUG writeSPIDER: File %s\n", filename.c_str());
#endif
    //locking
    struct flock fl;

    fl.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
    fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
    fl.l_start  = 0;        /* Offset from l_whence         */
    fl.l_len    = 0;        /* length, 0 = to EOF           */
    fl.l_pid    = getpid(); /* our PID                      */


    /*
     * BLOCK HEADER IF NEEDED
     */
    fl.l_type   = F_WRLCK;
    fcntl(fileno(fimg), F_SETLKW, &fl); /* locked */
    if(mode==WRITE_OVERWRITE || mode==WRITE_APPEND)//header must change
        fwrite( header, offset, 1, fimg );

    char* fdata = (char *) askMemory(datasize);
    //think about writing in several chucks

    //write only once, ignore select_img
    if ( NSIZE(data) == 1 && mode==WRITE_OVERWRITE)
    {
    	castPage2Datatype(MULTIDIM_ARRAY(data), fdata, Float, datasize_n);
        fwrite( fdata, datasize, 1, fimg );
    }

    else
    {
        if(mode==WRITE_APPEND)
            fseek( fimg, 0, SEEK_END);
        else if(mode==WRITE_REPLACE)
            fseek( fimg,offset + (offset+datasize)*select_img, SEEK_SET);

        // SJORS 30Oct12: I am completely unsure whether the code below will actually work....
        // Let's just raise an error an go out...
        REPORT_ERROR("writeSPIDER append/replace writing of SPIDER stacks not implemented yet....");
    }
    //I guess I do not need to unlock since we are going to close the file
    fl.l_type   = F_UNLCK;
    fcntl(fileno(fimg), F_SETLK, &fl); /* unlocked */

    freeMemory(fdata, datasize);
    freeMemory(header, (int)labbyt*sizeof(char));

    return(0);
}
#endif

