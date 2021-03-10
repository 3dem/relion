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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/*
        Based on rwMRC.h
        Header file for reading and writing MRC files
        Format: 3D crystallographic image file format for the MRC package
        Author: Bernard Heymann
        Created: 19990321 Modified: 20030723
*/

#ifndef RWMRC_H
#define RWMRC_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define SAFESET(var, value) (var = (std::isfinite(value)) ? value : var)

#define MRCSIZE 1024 // Minimum size of the MRC header (when nsymbt = 0)

///@defgroup MRC MRC File format
///@ingroup ImageFormats

/** MRC Old Header
  * @ingroup MRC
*/
struct MRCheadold
{
	// file header for MRC data
	int nx;              //  0   0       image size
	int ny;              //  1   4
	int nz;              //  2   8
	int mode;            //  3           0=uchar,1=short,2=float
	int nxStart;         //  4           unit cell offset
	int nyStart;         //  5
	int nzStart;         //  6
	int mx;              //  7           unit cell size in voxels
	int my;              //  8
	int mz;              //  9
	float a;             // 10   40      cell dimensions in A
	float b;             // 11
	float c;             // 12
	float alpha;         // 13           cell angles in degrees
	float beta;          // 14
	float gamma;         // 15
	int mapc;            // 16           column axis
	int mapr;            // 17           row axis
	int maps;            // 18           section axis
	float amin;          // 19           minimum density value
	float amax;          // 20   80      maximum density value
	float amean;         // 21           average density value
	int ispg;            // 22           space group number
	int nsymbt;          // 23           bytes used for sym. ops. table
	float extra[29];     // 24           user-defined info
	float xOrigin;       // 53           phase origin in pixels
	float yOrigin;       // 54
	int nlabl;           // 55           number of labels used
	char labels[10][80]; // 56-255       10 80-character labels
};

/** MRC Header
  * @ingroup MRC
*/
struct MRChead
{
	// file header for MRC data
	int nx;              //  0   0       image size
	int ny;              //  1   4
	int nz;              //  2   8
	int mode;            //  3           0=char,1=short,2=float
	int nxStart;         //  4           unit cell offset
	int nyStart;         //  5
	int nzStart;         //  6
	int mx;              //  7           unit cell size in voxels
	int my;              //  8
	int mz;              //  9
	float a;             // 10   40      cell dimensions in A
	float b;             // 11
	float c;             // 12
	float alpha;         // 13           cell angles in degrees
	float beta;          // 14
	float gamma;         // 15
	int mapc;            // 16           column axis
	int mapr;            // 17           row axis
	int maps;            // 18           section axis
	float amin;          // 19           minimum density value
	float amax;          // 20   80      maximum density value
	float amean;         // 21           average density value
	int ispg;            // 22           space group number
	int nsymbt;          // 23           bytes used for sym. ops. table
	float extra[25];     // 24           user-defined info
	float xOrigin;       // 49           phase origin in pixels
	float yOrigin;       // 50
	float zOrigin;       // 51
	char map[4];         // 52       identifier for map file ("MAP ")
	char machst[4];      // 53           machine stamp
	float arms;          // 54       RMS deviation
	int nlabl;           // 55           number of labels used
	char labels[800];    // 56-255       10 80-character labels
} ;

// For determination of machine stamp
#define UNKNOWN_SYSTEM 0
#define BIGIEEE 1
#define LITTLEIEEE 2
#define LITTLEVAX 3
int systype()
{
	char *test = (char*)askMemory(12);
	int *itest = (int*)test;
	float *ftest = (float*)test;
	memcpy(test, "jbh     ", 8);

	int type = UNKNOWN_SYSTEM;
	if (itest[0] == 1784834080 && fabs(ftest[0] - 6.84272e+25) < 1e+21)
		type = BIGIEEE;
	if (itest[0] == 543711850 && fabs(ftest[0] - 1.96837e-19) < 1e-23)
		type = LITTLEIEEE;

	freeMemory(test, 12);

	return type;
}

// I/O prototypes
/** MRC Reader
  * @ingroup MRC
*/
int readMRC(long int img_select, bool isStack=false, const FileName &name="")
{
#undef DEBUG
//#define DEBUG
#ifdef DEBUG
	printf("DEBUG readMRC: Reading MRC file\n");
#endif

	MRChead *header = (MRChead*)askMemory(sizeof(MRChead));
	if (fread(header, MRCSIZE, 1, fimg) < 1)
		REPORT_ERROR("rwMRC: error in reading header of image " + name);

	// Determine byte order and swap bytes if from little-endian machine
	swap = 0;
	char *b = (char*) header;
	int i;
	if ((abs(header->mode) > SWAPTRIG) || (abs(header->nx) > SWAPTRIG))
	{
#ifdef DEBUG
		fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
#endif
		swap = 1;
		int extent = MRCSIZE - 800; // exclude labels from swapping
		for (i = 0; i < extent; i += 4)
			swapbytes(b+i, 4);
	}

	long int _xDim,_yDim,_zDim;
	long int _nDim;
	_xDim = (int) header->nx;
	_yDim = (int) header->ny;
	_zDim = (int) header->nz;
	_nDim = 1;

	if(isStack)
	{
		_nDim = (long int)_zDim;
		_zDim = 1;
		replaceNsize = _nDim;
		std::stringstream Num;
		std::stringstream Num2;
		if (img_select >= (int)_nDim) // img_select starts from 0, while _nDim from 1
		{
			Num << (img_select + 1);
			Num2 << _nDim;
			REPORT_ERROR((std::string)"readMRC: Image number " + Num.str() + " exceeds stack size " + Num2.str() + " of image " + name);
		}
	}
	else
	{
	    replaceNsize = 0;
	}

	// Map the parameters
	if (isStack && img_select == -1)
		_zDim = 1;
	else if (isStack && img_select!=-1)
		_zDim = _nDim = 1;
	else
		_nDim = 1;

	data.setDimensions(_xDim, _yDim, _zDim, _nDim);

	DataType datatype;

	switch (header->mode)
	{
	case 0:
		datatype = SChar; // Changed to SIGNED in RELION 3.1 to be compatible with the official speficifation and SerialEM 
		break;
	case 1:
		datatype = SShort;
		break;
	case 2:
		datatype = Float;
		break;
	case 3:
		REPORT_ERROR("readMRC: only real-space images may be read into RELION.");
	case 4:
		REPORT_ERROR("readMRC: only real-space images may be read into RELION.");
	case 6:
		datatype = UShort;
		break;
	case 12:
		datatype = Float16;
		break;
	case 101:
		// This is SerialEM's non-standard extension.
		// https://bio3d.colorado.edu/imod/doc/mrc_format.txt
		// http://bio3d.colorado.edu/SerialEM/hlp/html/hidd_k2_save_options.htm
		if (_xDim % 2 == 1 && _yDim % 2 == 1)
			REPORT_ERROR("Currently we support 4-bit MRC (mode 101) only when nx * ny is an even number.");
		datatype = UHalf;
		break;
	default:
		REPORT_ERROR((std::string)"readMRC: unsupported MRC mode " + integerToString(header->mode));
	}
	offset = MRCSIZE + header->nsymbt;

	MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN,(RFLOAT)header->amin);
	MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX,(RFLOAT)header->amax);
	MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG,(RFLOAT)header->amean);
	MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV,(RFLOAT)header->arms);
	MDMainHeader.setValue(EMDL_IMAGE_DATATYPE,(int)datatype);

	if (header->mx && header->a != 0) //ux
		MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, (RFLOAT)header->a / header->mx);
	if (header->my && header->b != 0) //yx
		MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, (RFLOAT)header->b / header->my);
	if (header->mz && header->c != 0)//zx
		MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, (RFLOAT)header->c / header->mz);

	if (isStack && dataflag<0) // Don't read the individual header and the data if not necessary
	{
		freeMemory(header, sizeof(MRChead));
		return 0;
	}

//#define DEBUG
#ifdef DEBUG
	MDMainHeader.write(std::cerr);
	MD[0].write(std::cerr);
#endif

	freeMemory(header, sizeof(MRChead));

	return readData(fimg, img_select, datatype, 0);
}

/** MRC Writer
  * @ingroup MRC
*/
int writeMRC(long int img_select, bool isStack=false, const int mode=WRITE_OVERWRITE, const DataType datatype=Unknown_Type) /* TODO: add type */
{
	MRChead *header = (MRChead *) askMemory(sizeof(MRChead));

	// Map the parameters
	strncpy(header->map, "MAP ", 4);
	// Set CCP4 machine stamp
	switch (systype())
	{
	case BIGIEEE:
		header->machst[0] = header->machst[1] = 17;
		break;
	case LITTLEIEEE:
		header->machst[0] = 68;
		header->machst[1] = 65;
		break;
	case LITTLEVAX:
		header->machst[0] = 34;
		header->machst[1] = 65;
		break;
	case UNKNOWN_SYSTEM:
		REPORT_ERROR("Unkown system type in writeMRC machine stamp determination.");
	default:
		break;
	}

	// FIXME TO BE DONE WITH rwCCP4!!
	//set_CCP4_machine_stamp(header->machst);
	long int Xdim = XSIZE(data);
	long int Ydim = YSIZE(data);
	long int Zdim = ZSIZE(data);
	long int Ndim = NSIZE(data);
	long int imgStart = 0;
	long int imgEnd = Ndim;
	if (img_select != -1)
	{
		imgStart = img_select;
		imgEnd = img_select + 1;
	}
	if (mode == WRITE_APPEND || mode == WRITE_REPLACE)
	{
		imgStart = 0;
		imgEnd = 1;
	}
	header->nx = Xdim;
	header->ny = Ydim;
	if (isStack)
		header->nz = Ndim;
	else
		header->nz = Zdim;

	// Convert T to datatype
	DataType output_type;
	if ((datatype == Unknown_Type && (typeid(T) == typeid(RFLOAT) ||
	                                  typeid(T) == typeid(float) ||
	                                  typeid(T) == typeid(int)))
           || datatype == Float)
	{
		header->mode = 2;
		output_type = Float;
	}
	else if ((datatype == Unknown_Type && (typeid(T) == typeid(unsigned char) ||
	                                       typeid(T) == typeid(signed char)))
	        || datatype == SChar)
	{
		header->mode = 0;
		output_type = SChar;
	}
	else if ((datatype == Unknown_Type && typeid(T) == typeid(signed short))
	        || datatype == SShort)
	{
		header->mode = 1;
		output_type = SShort;
	}
	else if ((datatype == Unknown_Type && typeid(T) == typeid(unsigned short))
	        || datatype == UShort)
	{
		header->mode = 6;
		output_type = UShort;
	}
	else if (datatype == Float16)
	{
		header->mode = 12;
		output_type = Float16;
	}
	else
		REPORT_ERROR(std::string("writeMRC(): invalid output data type. datatype = ") + integerToString(datatype));

	//Set this to zero till we decide if we want to update it
	header->mx = header->nx; //(int) (ua/ux + 0.5);
	header->my = header->ny; //(int) (ub/uy + 0.5);
	header->mz = header->nz; //(int) (uc/uz + 0.5);
	header->mapc = 1;
	header->mapr = 2;
	header->maps = 3;
	RFLOAT aux,aux2;

	// TODO: fix this!
	header->a = header->nx; // ua;
	header->b = header->ny; // ub;
	header->c = header->nz; // uc;
	header->alpha = (float)90.;
	header->beta = (float)90.;
	header->gamma = (float)90.;
	header->xOrigin = (float)0.;
	header->yOrigin = (float)0.;
	header->zOrigin = (float)0.;
	header->nxStart = (int)0;
	header->nyStart = (int)0;
	header->nzStart = (int)0;

	if (!MDMainHeader.isEmpty())
	{
		if (MDMainHeader.getValue(EMDL_IMAGE_STATS_MIN, aux))
			header->amin = (float)aux;
		else
			header->amin = (float)data.computeMin();

		if (MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX, aux))
			header->amax = (float)aux;
		else
			header->amax = (float)data.computeMax();

		if (MDMainHeader.getValue(EMDL_IMAGE_STATS_AVG, aux))
			header->amean = (float)aux;
		else
			header->amean = (float)data.computeAvg();

		if (MDMainHeader.getValue(EMDL_IMAGE_STATS_STDDEV, aux))
			header->arms = (float)aux;
		else
			header->arms = (float)data.computeStddev();

		//if(MDMainHeader.getValue(EMDL_ORIENT_ORIGIN_X, aux))
		//	SAFESET(header->nxStart,(int)(aux-0.5));

		if (MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, aux2)) //header is init to zero
		{
			SAFESET(header->xOrigin, (float)(header->nxStart * aux2));
			SAFESET(header->a, (float)aux2 * header->nx);
		}

		//if (MDMainHeader.getValue(EMDL_ORIENT_ORIGIN_Y, aux))
		//	SAFESET(header->nyStart,(int)(aux-0.5));

		if (MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Y, aux2)) //header is init to zero
		{
			SAFESET(header->yOrigin, (float)(header->nyStart * aux2));
			SAFESET(header->b, (float)aux2 * header->ny);
		}

		//if (MDMainHeader.getValue(EMDL_ORIENT_ORIGIN_Z, aux))
		//	SAFESET(header->nzStart,(int)(aux-0.5));

		if (MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Z, aux2)) //header is init to zero
		{
			SAFESET(header->zOrigin, (float)(header->nzStart * aux2));
			SAFESET(header->c, (float)aux2 * header->nz);
		}

	}

	header->nsymbt = 0;

	//Create label "Relion version    date time"
#define MRC_LABEL_LEN 80
	header->nlabl = 1;

	char label[MRC_LABEL_LEN] = "Relion ";
	time_t rawtime;
	struct tm * timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

#ifdef PACKAGE_VERSION
	strcat(label,PACKAGE_VERSION);
#endif
	strcat(label, "   ");
	strftime(label + strlen(label), MRC_LABEL_LEN - strlen(label), "%d-%b-%y  %R:%S", timeinfo);
	strncpy(header->labels, label, MRC_LABEL_LEN);

	//strncpy(header->labels, p->label.c_str(), 799);

	offset = MRCSIZE + header->nsymbt;
	size_t datasize, datasize_n;
	datasize_n = Xdim * Ydim * Zdim;
	datasize = datasize_n * gettypesize(output_type);

//#define DEBUG
#ifdef DEBUG
	printf("DEBUG rwMRC: Offset = %ld,  Datasize_n = %ld\n", offset, datasize_n);
#endif

	// For multi-image files
	if (mode == WRITE_APPEND && isStack)
	{
		header->nz = replaceNsize + 1;
	}
	//else header-> is correct

	//locking
	struct flock fl;

	fl.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
	fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
	fl.l_start  = 0;        /* Offset from l_whence         */
	fl.l_len    = 0;        /* length, 0 = to EOF           */
	fl.l_pid    = getpid(); /* our PID                      */

	//BLOCK
	fl.l_type   = F_WRLCK;
	fcntl(fileno(fimg), F_SETLKW, &fl); /* locked */

	// Write header
	if(mode == WRITE_OVERWRITE || mode == WRITE_APPEND)
		fwrite(header, MRCSIZE, 1, fimg);
	freeMemory(header, sizeof(MRChead));

	//write only once, ignore select_img
	char* fdata = (char*)askMemory(datasize);
	//think about writing in several chunks

	if ( NSIZE(data) == 1 && mode==WRITE_OVERWRITE)
	{
		castPage2Datatype(MULTIDIM_ARRAY(data), fdata, output_type, datasize_n);
		fwrite(fdata, datasize, 1, fimg);
	}
	else
	{
		if (mode == WRITE_APPEND)
			fseek(fimg, 0, SEEK_END);
		else if (mode == WRITE_REPLACE)
		{
			fseek(fimg, offset + datasize * img_select, SEEK_SET);
		}

		for (size_t i = imgStart; i < imgEnd; i++)
		{
			castPage2Datatype(MULTIDIM_ARRAY(data) + i * datasize_n, fdata, output_type, datasize_n);
			fwrite(fdata, datasize, 1, fimg);
		}
	}

	// Unlock the file
	fl.l_type = F_UNLCK;
	fcntl(fileno(fimg), F_SETLK, &fl); /* unlocked */

	freeMemory(fdata, datasize);

	return(0);
}
#endif
