/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres", "Takanori Nakane"
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
/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef IMAGE_H
#define IMAGE_H

#include <cstdint>
#include <typeinfo>
#include <fcntl.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <tiffio.h>
#include "src/funcs.h"
#include "src/memory.h"
#include "src/filename.h"
#include "src/multidim_array.h"
#include "src/transformations.h"
#include "src/metadata_table.h"
#include "src/fftw.h"
#include "src/float16.h"

/// @defgroup Images Images
//@{

/** Data type.
 * This class defines the datatype of the data inside this image.
 */
typedef enum
{
	Unknown_Type = 0, // Undefined data type
	UChar = 1,        // Unsigned character or byte type
	SChar = 2,        // Signed character (for CCP4)
	UShort = 3,       // Unsigned integer (2-byte)
	SShort = 4,        // Signed integer (2-byte)
	UInt = 5,         // Unsigned integer (4-byte)
	Int = 6,          // Signed integer (4-byte)
	Long = 7,         // Signed integer (4 or 8 byte, depending on system)
	Float = 8,        // Floating point (4-byte)
	Double = 9,       // Double precision floating point (8-byte)
	Boolean = 10,     // Boolean (1-byte?)
	UHalf = 11,       // Signed 4-bit integer (SerialEM extension)
	Float16 = 12,     // Half precision floating point (2-byte)
	LastEntry = 15    // This must be the last entry
} DataType;

/** Write mode
 * This class defines the writing behavior.
 */
typedef enum
{
	WRITE_OVERWRITE, //forget about the old file and overwrite it
	WRITE_APPEND,	 //append and object at the end of a stack, so far can not append stacks
	WRITE_REPLACE,	 //replace a particular object by another
	WRITE_READONLY	 //only can read the file
} WriteMode;

extern "C" {
	typedef struct TiffInMemory
	{
		unsigned char *buf;
		tsize_t size;
		toff_t pos;
	} TiffInMemory;

	static tsize_t TiffInMemoryReadProc(thandle_t handle, tdata_t buf, tsize_t read_size)
	{
		TiffInMemory *tiff_handle = (TiffInMemory*)handle;
#ifdef TIFF_DEBUG
		std::cout << "TiffInMemoryReadProc: read_size = " << read_size << " cur_pos = " << tiff_handle->pos << " buf_size = " << tiff_handle->size << std::endl;
#endif
		if (tiff_handle->pos + read_size >= tiff_handle->size)
			REPORT_ERROR("TiffInMemoryReadProc: seeking beyond the end of the buffer.");

		memcpy(buf, tiff_handle->buf + tiff_handle->pos, read_size);
		tiff_handle->pos += read_size;

		return read_size;
	}

	static tsize_t TiffInMemoryWriteProc(thandle_t handle, tdata_t buf, tsize_t write_size)
	{
#ifdef TIFF_DEBUG
		REPORT_ERROR("TiffInMemoryWriteProc: Not implemented.");
#endif

		return -1;
	}

	static toff_t TiffInMemorySeekProc(thandle_t handle, toff_t offset, int whence)
	{
		TiffInMemory *tiff_handle = (TiffInMemory*)handle;
#ifdef TIFF_DEBUG
		std::cout << "TiffInMemorySeekProc: offset = " << offset << " cur_pos = " << tiff_handle->pos << " buf_size = " << tiff_handle->size << std::endl;
#endif
		switch (whence)
		{
			case SEEK_SET:
				tiff_handle->pos = 0;
				break;
			case SEEK_CUR:
				tiff_handle->pos += offset;
				break;
			case SEEK_END:
				REPORT_ERROR("TIFFInMemorySeekProc: SEEK_END is not supported.");
				// break; // intentional to suppress compiler warnings.
		}

		if (tiff_handle->pos >= tiff_handle->size)
			REPORT_ERROR("TIFFInMemorySeekProc: seeking beyond the end of the buffer.");

		return 0;
	}

	static int TiffInMemoryCloseProc(thandle_t handle)
	{
#ifdef TIFF_DEBUG
		std::cout << "TiffInMemoryCloseProc" << std::endl;
#endif
		return 0;
	}

	static toff_t TiffInMemorySizeProc(thandle_t handle)
	{
#ifdef TIFF_DEBUG
		std::cout << "TiffInMemorySizeProc" << std::endl;
#endif
		return ((TiffInMemory*)handle)->size;
	}

	static int TiffInMemoryMapFileProc(thandle_t handle, tdata_t *base, toff_t *size)
	{
		TiffInMemory *tiff_handle = (TiffInMemory*)handle;
#ifdef TIFF_DEBUG
		std::cout << "TiffInMemoryMapFileProc" << std::endl;
#endif

		*base = tiff_handle->buf;
		*size = tiff_handle->size;

		return 1;
	}

	static void TiffInMemoryUnmapFileProc(thandle_t handle, tdata_t base, toff_t size)
 	{
#ifdef TIFF_DEBUG
		std::cout << "TiffInMemoryUnmapFileProc" << std::endl;
#endif

		return;
	}
}

/** File handler class
 * This struct is used to share the File handlers with Image Collection class
 */
class fImageHandler
{
public:
	FILE*	  fimg;	// Image File handler
	FILE*	  fhed;	// Image File header handler
	TIFF*	  ftiff;
	FileName  ext_name; // Filename extension
	bool	  exist;    // Shows if the file exists
	bool	  isTiff;   // Shows if this is a TIFF file

	/** Empty constructor
	 */
	fImageHandler()
	{
		fimg=NULL;
		fhed=NULL;
		ftiff=NULL;
		ext_name="";
		exist=false;
		isTiff=false;
	}

	/** Destructor: closes file (if it still open)
	 */
	~fImageHandler()
	{
		closeFile();
	}

	void openFile(const FileName &name, int mode = WRITE_READONLY)
	{
		// This is thread-safe in C++11.
		// https://stackoverflow.com/questions/8102125/is-local-static-variable-initialization-thread-safe-in-c11
		static const size_t bufferSize = []() -> size_t {
			char * bufferString = getenv("RELION_STACK_BUFFER");
			if (bufferString != NULL)
				return std::atoi(bufferString);
    			else
				return SIZE_MAX;
		}();

		// Close any file that was left open in this handler
		if (!(fimg ==NULL && fhed == NULL))
			closeFile();

		FileName fileName, headName = "";
		// get the format, checking for possible format specifier before suffix
		// getFileFormat("file.spi") will return "spi"
		// getFileFormat("file.spi:mrc") will return "mrc"
		// getFileFormat("file") will return ""
		ext_name = name.getFileFormat();

		long int dump;
		name.decompose(dump, fileName);
		// Subtract 1 to have numbering 0...N-1 instead of 1...N
		if (dump > 0)
			dump--;

		// create the filename from a possible input format specifier (file.spi:mrc means "it's called .spi, but it's really a .mrc")
		// file.spi:mrc -> file.spi
		fileName = fileName.removeFileFormat();

		size_t found = fileName.find_first_of("%");
		if (found!=std::string::npos)
		  fileName = fileName.substr(0, found) ;

		exist = exists(fileName);

		std::string wmChar;

		switch (mode)
		{
		case WRITE_READONLY:
			if (!exist)
				REPORT_ERROR((std::string) "Cannot read file " + fileName + " It does not exist" );
			wmChar = "r";
			break;
		case WRITE_OVERWRITE:
			wmChar = "w";
			break;
		case WRITE_APPEND:
			if (exist)
				wmChar = "r+";
			else
				wmChar = "w+";
			break;
		case WRITE_REPLACE:
			wmChar = "r+";
			break;
		}

		if (ext_name.contains("img") || ext_name.contains("hed"))
		{
			fileName = fileName.withoutExtension();
			headName = fileName.addExtension("hed");
			fileName = fileName.addExtension("img");
		}
		else if(ext_name=="")
		{
			ext_name="spi"; // SPIDER is default format if none is specified
			fileName = fileName.addExtension(ext_name);
		}

		isTiff = ext_name.contains("tif");

		// Open image file
		if (isTiff) 
		{
			if (mode != WRITE_READONLY)
				REPORT_ERROR((std::string)"TIFF is supported only for reading");

			if ((ftiff = TIFFOpen(fileName.c_str(), "r")) == NULL)
				REPORT_ERROR((std::string)"Image::openFile cannot open: " + name);
		}
		else
		{
			if ((fimg  = fopen(fileName.c_str(), wmChar.c_str())) == NULL)
				REPORT_ERROR((std::string)"Image::openFile cannot open: " + name);

			if (ext_name=="mrcs" && wmChar=="r")
			{
				if (bufferSize < SIZE_MAX)
				{
					if (bufferSize == 0)
					{
						//disabling buffered IO for mrcs stacks to improve random IO behavior
						setvbuf(fimg, NULL, _IONBF, 0);
					}
					else
					{
						//set custom buffer size
						setvbuf(fimg, NULL, _IOFBF, bufferSize);
					}
				}
			}
		}

		if (headName != "")
		{
			if ((fhed = fopen(headName.c_str(), wmChar.c_str())) == NULL)
				REPORT_ERROR((std::string)"Image::openFile cannot open: " + headName);
		}
		else
			fhed = NULL;

	}

	void closeFile()
	{
		ext_name="";
		exist=false;

		// Check whether the file was closed already
		if (fimg == NULL && fhed == NULL && ftiff == NULL)
			return;

		if (isTiff && ftiff != NULL) {
			TIFFClose(ftiff);
			ftiff = NULL;
		}

		if (!isTiff && fclose(fimg) != 0)
			REPORT_ERROR((std::string)"Can not close image file ");
		else
			fimg = NULL;

		if (fhed != NULL &&  fclose(fhed) != 0)
			REPORT_ERROR((std::string)"Can not close header file ");
		else
			fhed = NULL;
	}

};

/** Returns memory size of datatype
 */
unsigned long gettypesize(DataType type);

/** Convert datatype string to datatypr enun */
int datatypeString2Int(std::string s);

/** Swapping trigger.
 * Threshold file z size above which bytes are swapped.
 */
#define SWAPTRIG	 65535

/** Template class for images.
 * The image class is the general image handling class.
 */
template<typename T>
class Image
{
public:
	MultidimArray<T> data; // The image data array
	MetaDataTable MDMainHeader; // metadata for the file

private:
	FileName filename; // File name
	FILE* fimg; // Image File handler
	FILE* fhed; // Image File header handler
	bool stayOpen; // To maintain the image file open after read/write
	int dataflag; // Flag to force reading of the data
	unsigned long i; // Current image number (may be > NSIZE)
	unsigned long offset; // Data offset
	int swap; // Perform byte swapping upon reading
	long int replaceNsize; // Stack size in the replace case
	bool _exists;  // does target file exists?
	// equal 0 is not exists or not a stack
	bool mmapOn; // Mapping when loading from file
	int mFd; // Handle the file in reading method and mmap
	size_t mappedSize; // Size of the mapped file

public:
	/** Empty constructor
	 *
	 * An empty image is created.
	 *
	 * @code
	 * Image<RFLOAT> I;
	 * @endcode
	 */
	Image()
	{
		mmapOn = false;
		clear();
		MDMainHeader.addObject();
	}

	/** Constructor with size
	 *
	 * A blank image (0.0 filled) is created with the given size. Pay attention
	 * to the dimension order: Y and then X.
	 *
	 * @code
	 * Image I(64,64);
	 * @endcode
	 */
	Image(long int Xdim, long int Ydim, long int Zdim=1, long int Ndim=1)
	{
		mmapOn = false;
		clear();
		data.resize(Ndim, Zdim, Ydim, Xdim);
		MDMainHeader.addObject();
	}

	/** Clear.
	 * Initialize everything to 0
	 */
	void clear()
	{
		if (mmapOn)
		{
			munmap(data.data-offset,mappedSize);
			close(mFd);
			data.data = NULL;
		}
		else
			data.clear();

		dataflag = -1;
		i = 0;
		filename = "";
		offset = 0;
		swap = 0;
		clearHeader();
		replaceNsize=0;
		mmapOn = false;
	}

	/** Clear the header of the image
	 */
	void clearHeader()
	{
		MDMainHeader.clear();
	}

	/** Destructor.
	 */
	~Image()
	{
		clear();
	}


	/** Specific read functions for different file formats
	  */
#include "src/rwSPIDER.h"
#include "src/rwMRC.h"
#include "src/rwIMAGIC.h"
#include "src/rwTIFF.h"

	/** Is this file an image
	 *
	 *	Check whether a real-space image can be read
	 *
	 */
	bool isImage(const FileName &name)
	{
		return !read(name, false);
	}

	/** Rename the image
	  */
	void rename (const FileName &name)
	{
		filename = name;
	}

	/** General read function
	 * you can read a single image from a single image file
	 * or a single image file from an stack, in the second case
	 * the select slide may come in the image name or in the select_img parameter
	 * file name takes precedence over select_img
	 * If -1 is given the whole object is read
	 * The number before @ in the filename is 1-indexed, while select_img is 0-indexed.
	 */
	int read(const FileName &name, bool readdata=true, long int select_img=-1, bool mapData = false, bool is_2D = false)
	{

		if (name == "")
			REPORT_ERROR("ERROR: trying to read image with empty file name!");
		int err = 0;
		fImageHandler hFile;
		hFile.openFile(name);
		err = _read(name, hFile, readdata, select_img, mapData, is_2D);
		// the destructor of fImageHandler will close the file

		// Negative errors are bad
		return err;
	}

	/** Read function from a file that has already been opened
	 *
	 */
	int readFromOpenFile(const FileName &name, fImageHandler &hFile, long int select_img, bool is_2D = false)
	{
		int err = 0;
		err = _read(name, hFile, true, select_img, false, is_2D);
		// Reposition file pointer for a next read
		rewind(fimg);
		return err;
	}

	/** General write function
	 * select_img= which slice should I replace
	 * overwrite = 0, append slice
	 * overwrite = 1 overwrite slice
	 *
	 * NOTE:
	 *	select_img has higher priority than the number before "@" in the name.
	 *	select_img counts from 0, while the number before "@" in the name from 1!
	 */
	void write(FileName name="",
	           long int select_img=-1,
	           bool isStack=false,
	           const WriteMode mode=WRITE_OVERWRITE,
	           const DataType datatype=Unknown_Type)
	{

		const FileName &fname = (name == "") ? filename : name;
		fImageHandler hFile;
		hFile.openFile(name, mode);
		_write(fname, hFile, select_img, isStack, mode, datatype);
		// the destructor of fImageHandler will close the file

	}

	/** Cast a page of data from type dataType to type Tdest
	 *	  input pointer  char *
	 */
	void castPage2T(char *page, T *ptrDest, DataType datatype, size_t pageSize )
	{
		switch (datatype)
		{
		case Unknown_Type:
			REPORT_ERROR("ERROR: datatype is Unknown_Type");
		case UChar:
			{
				if (typeid(T) == typeid(unsigned char))
					memcpy(ptrDest, page, pageSize * sizeof(T));
				else
				{
					unsigned char *ptr = (unsigned char *)page;
					for (size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case SChar:
			{
				if (typeid(T) == typeid(signed char))
				{
					memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					signed char *ptr = (signed char *)page;
					for (size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case UShort:
			{
				if (typeid(T) == typeid(unsigned short))
				{
					memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					unsigned short *ptr = (unsigned short *)page;
					for(size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case SShort:
			{
				if (typeid(T) == typeid(short))
				{
					memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					short *ptr = (short *)page;
					for(size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case UInt:
			{
				if (typeid(T) == typeid(unsigned int))
				{
					memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					unsigned int *ptr = (unsigned int *)page;
					for(size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case Int:
			{
				if (typeid(T) == typeid(int))
				{
					memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					int *ptr = (int *)page;
					for(size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case Long:
			{
				if (typeid(T) == typeid(long))
				{
					memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					long *ptr = (long *)page;
					for(size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case Float:
			{
				if (typeid(T) == typeid(float))
				{
				memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					float *ptr = (float *)page;
					for(size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case Double:
			{
				if (typeid(T) == typeid(RFLOAT))
				{
					memcpy(ptrDest, page, pageSize * sizeof(T));
				}
				else
				{
					RFLOAT *ptr = (RFLOAT *)page;
					for(size_t i = 0; i < pageSize; i++)
						ptrDest[i] = (T)ptr[i];
				}
				break;
			}
		case Float16:
			{
				float16 *ptr = (float16 *)page;
				for(size_t i = 0; i < pageSize; i++)
					ptrDest[i] = (T)half2float(ptr[i]);
				break;
			}
		case UHalf:
			{
				if (pageSize % 2 != 0) REPORT_ERROR("Logic error in castPage2T; for UHalf, pageSize must be even.");

				for(size_t i = 0, ilim = pageSize / 2; i < ilim; i++)
				{
					// Here we are assuming the fill-order is LSB2MSB according to IMOD's
					// iiProcessReadLine() in libiimod/mrcsec.c.
					// The default fill-order in the TIFF specification is MSB2LSB
					// but IMOD assumes LSB2MSB even for TIFF.
					// See IMOD's iiTIFFCheck() in libiimod/iitif.c.
					ptrDest[i * 2 ] = (T)(page[i] & 15); // 1111 = 1+2+4+8 = 15
					ptrDest[i * 2 + 1] = (T)((page[i] >> 4) & 15);
				}
				break;
			}
		default:
				{
					std::cerr<<"Datatype= "<<datatype<<std::endl;
					REPORT_ERROR(" ERROR: cannot cast datatype to T");
				}
			}

	}

	/** Cast page from T to datatype
	 *  input pointer char *
	 */
	void castPage2Datatype(T *srcPtr, char *page, DataType datatype, size_t pageSize)
	{
		switch (datatype)
		{
		case Float:
			{
				if (typeid(T) == typeid(float))
				{
					memcpy(page, srcPtr, pageSize*sizeof(T));
				}
				else
				{
					float *ptr = (float *)page;
					for (size_t i = 0; i < pageSize; i++)
						ptr[i] = (float)srcPtr[i];
				}
				break;
			}
		case Double:
			{
				if (typeid(T) == typeid(RFLOAT))
				{
					memcpy(page, srcPtr, pageSize*sizeof(T));
				}
				else
				{
					RFLOAT *ptr = (RFLOAT *)page;
					for (size_t i = 0; i < pageSize; i++)
						ptr[i] = (RFLOAT)srcPtr[i];
				}
				break;
			}
		case Float16:
			{
				float16 *ptr = (float16 *)page;
				for (size_t i = 0; i < pageSize; i++)
					ptr[i] = float2half((float)srcPtr[i]);
				break;
			}
		case SShort: 
			{
				if (typeid(T) == typeid(short))
				{
					memcpy(page, srcPtr, pageSize*sizeof(T));
				}
				else
				{
					short *ptr = (short *)page;
					for (size_t i = 0; i < pageSize; i++)
						ptr[i] = (short)srcPtr[i];
				}
				break;
			}
		case UShort:
			{
				if (typeid(T) == typeid(unsigned short))
				{
					memcpy(page, srcPtr, pageSize*sizeof(T));
				}
				else
				{
					unsigned short *ptr = (unsigned short *)page;
					for (size_t i = 0; i < pageSize; i++)
						ptr[i] = (unsigned short)srcPtr[i];
				}
				break;
			}
		case UChar:
			{
				if (typeid(T) == typeid(unsigned char))
				{
					memcpy(page, srcPtr, pageSize*sizeof(T));
				}
				else
				{
					unsigned char *ptr = (unsigned char *)page;
					for (size_t i = 0; i < pageSize; i++)
						ptr[i] = (unsigned char)srcPtr[i];
				}
				break;
			}
		default:
				{
					std::cerr << "outputDatatype= " << datatype << std::endl;
					REPORT_ERROR(" ERROR: cannot cast T to outputDatatype");
				}
			}
	}

	/** Check file Datatype is same as T type to use mmap.
	 */
	bool checkMmapT(DataType datatype)
	{

		switch (datatype)
		{
		case Unknown_Type:
			REPORT_ERROR("ERROR: datatype is Unknown_Type");
		case UChar:
			{
				if (typeid(T) == typeid(unsigned char))
					return 1;
				else
					return 0;
			}
		case SChar:
			{
				if (typeid(T) == typeid(signed char))
					return 1;
				else
					return 0;
			}
		case UShort:
			{
				if (typeid(T) == typeid(unsigned short))
					return 1;
				else
					return 0;
			}
		case SShort:
			{
				if (typeid(T) == typeid(short))
					return 1;
				else
					return 0;
			}
		case UInt:
			{
				if (typeid(T) == typeid(unsigned int))
					return 1;
				else
					return 0;
			}
		case Int:
			{
				if (typeid(T) == typeid(int))
					return 1;
				else
					return 0;
			}
		case Long:
			{
				if (typeid(T) == typeid(long))
					return 1;
				else
					return 0;
			}
		case Float:
			{
				if (typeid(T) == typeid(float))
					return 1;
				else
					return 0;
			}
		case Double:
			{
				if (typeid(T) == typeid(RFLOAT)) // TODO: CHECKME: Is this safe with single precision build?
					return 1;
				else
					return 0;
			}
		default:
			{
				std::cerr << "Datatype= " << datatype << std::endl;
				REPORT_ERROR(" ERROR: cannot cast datatype to T");
			}
		}
		//  int * iTemp = (int*) map;
		//  ptrDest = reinterpret_cast<T*> (iTemp);
	}

	/** Write an entire page as datatype
	 *
	 * A page of datasize_n elements T is cast to datatype and written to fimg
	 * The memory for the casted page is allocated and freed internally.
	 */
	void writePageAsDatatype(FILE * fimg, DataType datatype, size_t datasize_n )
	{
		size_t datasize = datasize_n * gettypesize(datatype);
		char * fdata = (char *) askMemory(datasize);
		castPage2Datatype(MULTIDIM_ARRAY(data), fdata, datatype, datasize_n);
		fwrite( fdata, datasize, 1, fimg );
		freeMemory(fdata, datasize);
	}

	/** Swap an entire page
	  * input pointer char *
	  */
	void swapPage(char * page, size_t pageNrElements, DataType datatype)
	{
		unsigned long datatypesize = gettypesize(datatype);
#ifdef DEBUG

		std::cerr<<"DEBUG swapPage: Swapping image data with swap= "
		<< swap<<" datatypesize= "<<datatypesize
		<< " pageNrElements " << pageNrElements
		<< " datatype " << datatype
		<<std::endl;
		;
#endif

		// Swap bytes if required
		if ( swap == 1 )
		{
			for (unsigned long i=0; i<pageNrElements; i+=datatypesize)
				swapbytes(page+i, datatypesize);
		}
		else if ( swap > 1 )
		{
			for (unsigned long i=0; i<pageNrElements; i+=swap)
				swapbytes(page+i, swap);
		}
	}

	/** Read the raw data
	  */
	int readData(FILE* fimg, long int select_img, DataType datatype, unsigned long pad)
	{
//#define DEBUG
#ifdef DEBUG
		std::cerr<<"entering readdata"<<std::endl;
		std::cerr<<" readData flag= "<<dataflag<<std::endl;
#endif

		if ( dataflag < 1 )
			return 0;

		size_t myoffset, readsize, readsize_n, pagemax = 1073741824; // 1 GB
		size_t datatypesize; // bytes
		size_t pagesize; // bytes
		if (datatype == UHalf)
		{
			if (YXSIZE(data) % 2 != 0) REPORT_ERROR("For UHalf, YXSIZE(data) must be even.");
			pagesize = ZYXSIZE(data) / 2;
		}
		else
		{
			datatypesize = gettypesize(datatype);
			pagesize = ZYXSIZE(data)*datatypesize;
		}
		size_t haveread_n=0; // number of pixels (not necessarily bytes!) processed so far

		//Multidimarray mmapOn is priority over image mmapOn
		if(data.mmapOn)
			mmapOn = false;

		if (datatype == UHalf) mmapOn = false;

		// Flag to know that data is not going to be mapped although mmapOn is true
		if (mmapOn && !checkMmapT(datatype))
		{
			std::cout << "WARNING: Image Class. File datatype and image declaration not compatible with mmap. Loading into memory." <<std::endl;
			mmapOn = false;
			mFd = -1;
		}

		if (mmapOn)
		{
			if ( NSIZE(data) > 1 )
			{
				REPORT_ERROR("Image Class::ReadData: mmap with multiple \
							 images file not compatible. Try selecting a unique image.");
			}

			fclose(fimg);

			//if ( ( mFd = open(filename.c_str(), O_RDWR, S_IREAD | S_IWRITE) ) == -1 )
			if ( ( mFd = open(filename.c_str(), O_RDWR, S_IRUSR | S_IWUSR) ) == -1 )
				REPORT_ERROR("Image Class::ReadData: Error opening the image file.");

			char * map;
			mappedSize = pagesize+offset;

			if ( (map = (char*) mmap(0,mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0)) == (void*) -1 )
				REPORT_ERROR("Image Class::ReadData: mmap of image file failed.");
			data.data = reinterpret_cast<T*> (map+offset);
		}
		else
		{
			// Reset select to get the correct offset
			if ( select_img < 0 )
				select_img = 0;

			char* page = NULL;

			// Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
			// if memory already allocated use it (no resize allowed)
			data.coreAllocateReuse();
			myoffset = offset + select_img*(pagesize + pad);
			//#define DEBUG

#ifdef DEBUG

			data.printShape();
			printf("DEBUG: Page size: %ld offset= %d \n", pagesize, offset);
			printf("DEBUG: Swap = %d  Pad = %ld  Offset = %ld\n", swap, pad, offset);
			printf("DEBUG: myoffset = %d select_img= %d \n", myoffset, select_img);
#endif

			if (pagesize > pagemax)
				page = (char *) askMemory(pagemax*sizeof(char));
			else
				page = (char *) askMemory(pagesize*sizeof(char));

			// Because we requested XYSIZE to be even for UHalf, this is always safe.
			int error_fseek = fseek(fimg, myoffset, SEEK_SET);
			if (error_fseek != 0)
				return -1;

			for (size_t myn=0; myn<NSIZE(data); myn++)
			{
				for (size_t myj=0; myj<pagesize; myj+=pagemax) //pagesize size of object
				{
					// Read next page. Divide pages larger than pagemax
					readsize = pagesize - myj;
					if ( readsize > pagemax )
						readsize = pagemax;

					if (datatype == UHalf)
						readsize_n = readsize * 2;
					else
						readsize_n = readsize/datatypesize;

#ifdef DEBUG
					std::cout << "NX = " << XSIZE(data) << " NY = " << YSIZE(data) << " NZ = " << ZSIZE(data) << std::endl;
					std::cout << "pagemax = " << pagemax << " pagesize = " << pagesize  << " readsize = " << readsize << " readsize_n = " << readsize_n << std::endl;
#endif

					//Read page from disc
					size_t result = fread( page, readsize, 1, fimg );
					if (result != 1)
						return -2;

					//swap per page
					if (swap)
						swapPage(page, readsize, datatype);
					// cast to T per page
					castPage2T(page, MULTIDIM_ARRAY(data) + haveread_n, datatype, readsize_n);
					haveread_n += readsize_n;
				}
				if ( pad > 0 )
				{
					//fread( padpage, pad, 1, fimg);
					error_fseek = fseek( fimg, pad, SEEK_CUR );
					if (error_fseek != 0)
						return -1;
				}
			}
			//if ( pad > 0 )
			//	  freeMemory(padpage, pad*sizeof(char));
			if ( page != NULL )
				freeMemory(page, pagesize*sizeof(char));

#ifdef DEBUG

			printf("DEBUG img_read_data: Finished reading and converting data\n");
#endif

		}
		return 0;
	}

	/** Data access
	 *
	 * This operator can be used to access the data multidimarray.
	 * In this way we could resize an image just by
	 * resizing its associated matrix or we could add two images by adding their
	 * matrices.
	 * @code
	 * I().resize(128, 128);
	 * I2() = I1() + I2();
	 * @endcode
	 */
	MultidimArray<T>& operator()()
	{
		return data;
	}
	const MultidimArray<T>& operator()() const
	{
		return data;
	}

	/** Pixel access
	*
	* This operator is used to access a pixel within a 2D image. This is a
	* logical access, so you could access to negative positions if the image
	* has been defined so (see the general explanation for the class).
	*
	* @code
	* std::cout << "Grey level of pixel (-3,-3) of the image = " << I(-3, -3)
	* << std::endl;
	*
	* I(-3, -3) = I(-3, -2);
	* @endcode
	*/
	T& operator()(int i, int j) const
	{
		return A2D_ELEM(data, i, j);
	}
	/** Set pixel
	 * (direct access) needed by swig
	 */
	void setPixel(int i, int j, T v)
	{
		IMGPIXEL(*this,i,j)=v;
	}

	/** Get pixel
	 * (direct acces) needed by swig
	 */
	T getPixel(int i, int j) const
	{
		return IMGPIXEL(*this,i,j);
	}

	/** Voxel access
	 *
	 * This operator is used to access a voxel within a 3D image. This is a
	 * logical access, so you could access to negative positions if the image
	 * has been defined so (see the general explanation for the class).
	 *
	 * @code
	 * std::cout << "Grey level of pixel (-3,-3, 1) of the volume = " << I(-3, -3, 1)
	 * << std::endl;
	 *
	 * I(-3, -3, 1) = I(-3, -2, 0);
	 * @endcode
	 */
	T& operator()(int k, int i, int j) const
	{
		return A3D_ELEM(data, k, i, j);
	}

	/** Get file name
	 *
	 * @code
	 * std::cout << "Image name = " << I.name() << std::endl;
	 * @endcode
	 */
	const FileName & name() const
	{
		return filename;
	}

	/** Get Image dimensions
	 */
	void getDimensions(int &Xdim, int &Ydim, int &Zdim, long int &Ndim) const
	{
		Xdim = XSIZE(data);
		Ydim = YSIZE(data);
		Zdim = ZSIZE(data);
		Ndim = NSIZE(data);
	}

	long unsigned int getSize() const
	{
		return NZYXSIZE(data);
	}

	/* Is there label in the main header */
	bool mainContainsLabel(EMDLabel label) const
	{
		return MDMainHeader.containsLabel(label);
	}

	/** Data type
		*
		* @code
		* std::cout << "datatype= " << dataType() << std::endl;
		* @endcode
		*/
	int dataType() const
	{
		int dummy;
		MDMainHeader.getValue(EMDL_IMAGE_DATATYPE, dummy);
		return dummy;
	}

	/** Sampling RateX
	*
	* @code
	* std::cout << "sampling= " << samplingRateX() << std::endl;
	* @endcode
	*/
	RFLOAT samplingRateX(const long int n = 0) const
	{
		RFLOAT dummy = 1.;
		MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, dummy);
		return dummy;
	}

	/** Sampling RateY
	*
	* @code
	* std::cout << "sampling= " << samplingRateY() << std::endl;
	* @endcode
	*/
	RFLOAT samplingRateY(const long int n = 0) const
	{
		RFLOAT dummy = 1.;
		MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Y, dummy);
		return dummy;
	}

	/** Set file name
	 */
	void setName(const FileName &_filename)
	{
		filename = _filename;
	}

	/* Set image statistics in the main header
	 *
	 */
	void setStatisticsInHeader()
	{
		RFLOAT avg,stddev,minval,maxval;
		data.computeStats(avg, stddev, minval, maxval);
		MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
		MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
		MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
		MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
	}

	void setSamplingRateInHeader(RFLOAT rate_x, RFLOAT rate_y = -1., RFLOAT rate_z = -1.)
	{
		MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, rate_x);
		if (rate_y < 0.)
			rate_y = rate_x;
		MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, rate_y);
		if (ZSIZE(data)>1)
		{
			if (rate_z < 0.)
				rate_z = rate_x;
			MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, rate_z);
		}
	}

	/** Show image properties
	  */
	friend std::ostream& operator<<(std::ostream& o, const Image<T>& I)
	{
		o << "Image type   : ";
			o << "Real-space image" << std::endl;

		o << "Reversed	   : ";
		if (I.swap)
			o << "TRUE"  << std::endl;
		else
			o << "FALSE" << std::endl;

		o << "Data type    : ";
		switch (I.dataType())
		{
		case Unknown_Type:
			o << "Undefined data type";
			break;
		case UChar:
			o << "Unsigned character or byte type";
			break;
		case SChar:
			o << "Signed character (for CCP4)";
			break;
		case UShort:
			o << "Unsigned integer (2-byte)";
			break;
		case SShort:
			o << "Signed integer (2-byte)";
			break;
		case UInt:
			o << "Unsigned integer (4-byte)";
			break;
		case Int:
			o << "Signed integer (4-byte)";
			break;
		case Long:
			o << "Signed integer (4 or 8 byte, depending on system)";
			break;
		case Float:
			o << "Floating point (4-byte)";
			break;
		case Double:
			o << "Double precision floating point (8-byte)";
			break;
		case Float16:
			o << "Half precision floating point (4-byte)";
			break;
		case Boolean:
			o << "Boolean (1-byte?)";
			break;
		case UHalf:
			o << "4-bit integer";
			break;
		}
		o << std::endl;

		o << "dimensions   : " << NSIZE(I()) << " x " << ZSIZE(I()) << " x " << YSIZE(I()) << " x " << XSIZE(I());
		o << "	(noObjects x slices x rows x columns)" << std::endl;
		return o;
	}

	/** Sum this object with other file and keep in this object
	  */
	void sumWithFile(const FileName &fn)
	{
		Image<T> aux;
		aux.read(fn);
		(*this)()+=aux();
	}

	int readTiffInMemory(void* buf, size_t size, bool readdata=true, long int select_img = -1,
	                     bool mapData = false, bool is_2D = false)
	{
		int err = 0;

		TiffInMemory handle;
		handle.buf = (unsigned char*)buf;
		handle.size = size;
		handle.pos = 0;
		// Check whether to read the data or only the header
		dataflag = ( readdata ) ? 1 : -1;

		// Check whether to map the data or not
		mmapOn = mapData;

		//Just clear the header before reading
		MDMainHeader.clear();
		MDMainHeader.addObject();

		TIFF* ftiff = TIFFClientOpen("in-memory-tiff", "r", (thandle_t)&handle,
		                             TiffInMemoryReadProc, TiffInMemoryWriteProc, TiffInMemorySeekProc,
		                             TiffInMemoryCloseProc, TiffInMemorySizeProc, TiffInMemoryMapFileProc,
		                             TiffInMemoryUnmapFileProc);
		err = readTIFF(ftiff, select_img, readdata, true, "in-memory-tiff");
		TIFFClose(ftiff);

		return err;
	}

private:
	int _read(const FileName &name, fImageHandler &hFile, bool readdata=true, long int select_img = -1,
			  bool mapData = false, bool is_2D = false)
	{
		int err = 0;

		// Check whether to read the data or only the header
		dataflag = ( readdata ) ? 1 : -1;

		// Check whether to map the data or not
		mmapOn = mapData;

		FileName ext_name = hFile.ext_name;
		fimg = hFile.fimg;
		fhed = hFile.fhed;

		long int dump;
		name.decompose(dump, filename);
		// Subtract 1 to have numbering 0...N-1 instead of 1...N
		if (dump > 0)
			dump--;
		filename = name;

		if (select_img == -1)
			select_img = dump;

#undef DEBUG
//#define DEBUG
#ifdef DEBUG
		std::cerr << "READ\n" <<
		"name="<<name <<std::endl;
		std::cerr << "ext= "<<ext_name <<std::endl;
		std::cerr << " now reading: "<< filename <<" dataflag= "<<dataflag
		<< " select_img "  << select_img << std::endl;
#endif
#undef DEBUG

		//Just clear the header before reading
		MDMainHeader.clear();
		MDMainHeader.addObject();

		if (ext_name.contains("spi") || ext_name.contains("xmp")  ||
			ext_name.contains("stk") || ext_name.contains("vol"))
			err = readSPIDER(select_img);
		else if (ext_name.contains("mrcs") || (is_2D && ext_name.contains("mrc")) || //mrc stack MUST go BEFORE plain MRC
				ext_name.contains("st")) //stk stack MUST go BEFORE plain st
			err = readMRC(select_img, true, name);
		else if (ext_name.contains("tif"))
			err = readTIFF(hFile.ftiff, select_img, readdata, true, name);
		else if (select_img >= 0 && ext_name.contains("mrc"))
			REPORT_ERROR("Image::read ERROR: stacks of images in MRC-format should have extension .mrcs; .mrc extensions are reserved for 3D maps.");
		else if (ext_name.contains("mrc")) // mrc 3D map
			err = readMRC(select_img, false, name);
		else if (ext_name.contains("img") || ext_name.contains("hed"))//
			err = readIMAGIC(select_img);//imagic is always an stack
		else if (ext_name.contains("dm"))
			REPORT_ERROR("The Digital Micrograph format (DM3, DM4) is not supported. You can convert it to MRC by other programs, for example, dm2mrc in IMOD.");
		else if (ext_name.contains("eer") || ext_name.contains("ecc"))
			REPORT_ERROR("BUG: EER movies should be handled by EERRenderer, not by Image.");
		else
			err = readSPIDER(select_img);

		// Negative errors are bad.
		return err;
	}

	void _write(const FileName &name, fImageHandler &hFile, long int select_img=-1,
	            bool isStack=false, const WriteMode mode=WRITE_OVERWRITE, const DataType datatype=Unknown_Type)
	{
		int err = 0;

		FileName ext_name = hFile.ext_name;
		fimg = hFile.fimg;
		fhed = hFile.fhed;
		_exists = hFile.exist;

		filename = name;

		long int aux;
		FileName filNamePlusExt(name);
		name.decompose(aux, filNamePlusExt);
		// Subtract 1 to have numbering 0...N-1 instead of 1...N
		if (aux > 0)
			aux--;

		if (select_img == -1)
			select_img = aux;

		size_t found = filNamePlusExt.find_first_of("%");

		std::string imParam = "";

		if (found!=std::string::npos)
		{
			imParam =  filNamePlusExt.substr(found+1).c_str();
			filNamePlusExt = filNamePlusExt.substr(0, found) ;
		}

		found = filNamePlusExt.find_first_of(":");
		if ( found!=std::string::npos)
			filNamePlusExt	 = filNamePlusExt.substr(0, found);

//#define DEBUG
#ifdef DEBUG
		std::cerr << "write" <<std::endl;
		std::cerr<<"extension for write= "<<ext_name<<std::endl;
		std::cerr<<"filename= "<<filename<<std::endl;
		std::cerr<<"mode= "<<mode<<std::endl;
		std::cerr<<"isStack= "<<isStack<<std::endl;
		std::cerr<<"select_img= "<<select_img<<std::endl;
#endif
#undef DEBUG
		// Check that image is not empty
		if (getSize() < 1)
			REPORT_ERROR("write Image ERROR: image is empty!");

		// CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???
		int Xdim, Ydim, Zdim;
		long int Ndim;
		this->getDimensions(Xdim,Ydim, Zdim, Ndim);

		Image<T> auxI;
		replaceNsize=0;//reset replaceNsize in case image is reused
		if(select_img == -1 && mode == WRITE_REPLACE)
			REPORT_ERROR("write: Please specify object to be replaced");
		else if(!_exists && mode == WRITE_REPLACE)
		{
			std:: stringstream replace_number;
			replace_number << select_img;
			REPORT_ERROR((std::string)"Cannot replace object number: "
						 + replace_number.str()
						 + " in file " +filename
						 + ". It does not exist");
		}
		else if (_exists && (mode == WRITE_REPLACE || mode == WRITE_APPEND))
		{
			auxI.dataflag = -2;
			auxI.read(filNamePlusExt,false);
			int _Xdim, _Ydim, _Zdim;
			long int _Ndim;
			auxI.getDimensions(_Xdim,_Ydim, _Zdim, _Ndim);
			replaceNsize=_Ndim;
			if(Xdim!=_Xdim ||
			   Ydim!=_Ydim ||
			   Zdim!=_Zdim
			  )
				REPORT_ERROR("write: target and source objects have different size");
			if(mode==WRITE_REPLACE && select_img>_Ndim)
				REPORT_ERROR("write: cannot replace image stack is not large enough");
			if(auxI.replaceNsize <1 &&
			   (mode==WRITE_REPLACE || mode==WRITE_APPEND))
				REPORT_ERROR("write: output file is not an stack");
		}
		else if(!_exists && mode==WRITE_APPEND)
		{
			;
		}
		else if (mode == WRITE_READONLY)//If new file we are in the WRITE_OVERWRITE mode
		{
			REPORT_ERROR( (std::string) "File " + name
						 + " opened in read-only mode. Cannot write.");
		}

		/*
		 * SELECT FORMAT
		 */
		if(ext_name.contains("spi") || ext_name.contains("xmp") ||
		   ext_name.contains("stk") || ext_name.contains("vol"))
			err = writeSPIDER(select_img, isStack, mode, datatype);
		else if (ext_name.contains("mrcs"))
			writeMRC(select_img, true, mode, datatype);
		else if (ext_name.contains("mrc"))
			writeMRC(select_img, false, mode, datatype);
		else if (ext_name.contains("img") || ext_name.contains("hed"))
			writeIMAGIC(select_img, mode);
		else
			err = writeSPIDER(select_img, isStack, mode, datatype);
		if ( err < 0 )
		{
			std::cerr << " Filename = " << filename << " Extension= " << ext_name << std::endl;
			REPORT_ERROR((std::string)"Error writing file "+ filename + " Extension= " + ext_name);
		}

		/* If initially the file did not exist, once the first image is written, then the file exists
		 */
		if (!_exists)
			hFile.exist = _exists = true;
	}
};

// Some image-specific operations

// For image normalisation
void normalise(Image<RFLOAT> &I,
               int bg_radius,
               RFLOAT white_dust_stddev,
               RFLOAT black_dust_stddev,
               bool do_ramp,
               bool is_helical_segment = false,
               RFLOAT helical_mask_tube_outer_radius_pix = -1.,
               RFLOAT tilt_deg = 0.,
               RFLOAT psi_deg = 0.);
void calculateBackgroundAvgStddev(Image<RFLOAT> &I,
                                  RFLOAT &avg,
                                  RFLOAT &stddev,
                                  int bg_radius,
                                  bool is_helical_segment = false,
                                  RFLOAT helical_mask_tube_outer_radius_pix = -1.,
                                  RFLOAT tilt_deg = 0.,
                                  RFLOAT psi_deg = 0.);
void subtractBackgroundRamp(Image<RFLOAT> &I,
                            int bg_radius,
                            bool is_helical_segment = false,
                            RFLOAT helical_mask_tube_outer_radius_pix = -1.,
                            RFLOAT tilt_deg = 0.,
                            RFLOAT psi_deg = 0.);

// For dust removal
void removeDust(Image<RFLOAT> &I, bool is_white, RFLOAT thresh, RFLOAT avg, RFLOAT stddev);

// for contrast inversion
void invert_contrast(Image<RFLOAT> &I);

// for image re-scaling
void rescale(Image<RFLOAT> &I, int mysize);

// for image re-windowing
void rewindow(Image<RFLOAT> &I, int mysize);

/// @defgroup ImageFormats Image Formats
/// @ingroup Images
// Functions belonging to this topic are commented in rw*.h
//@}

#define GREYSCALE 0
#define BLACKGREYREDSCALE 1
#define BLUEGREYWHITESCALE 2
#define BLUEGREYREDSCALE 3
#define RAINBOWSCALE 4
#define CYANBLACKYELLOWSCALE 5

void getImageContrast(MultidimArray<RFLOAT> &image, RFLOAT &minval, RFLOAT &maxval, RFLOAT &sigma_contrast);

inline void greyToRGB(const int color_scheme, const unsigned char grey, unsigned char &red, unsigned char &green, unsigned char &blue)
{
	switch (color_scheme)
	{
	case GREYSCALE:
		red = green = blue = grey;
		break;
	case BLACKGREYREDSCALE:
		if (grey >= 128) { red = 255; blue = green = FLOOR((RFLOAT)(255 - grey)*2); }
		else { red = green = blue = FLOOR((RFLOAT)(grey*2.)); }
		break;
	case BLUEGREYWHITESCALE:
		if (grey >= 128) { red = green = blue = FLOOR((RFLOAT)((grey - 128) * 2)); }
		else { red = 0; blue = green = FLOOR((RFLOAT)(255 - 2 * grey)); }
		break;
	case BLUEGREYREDSCALE:
	{
		const RFLOAT a = grey / 85.0; // group
		const int X = FLOOR(a);	//this is the integer part
		const unsigned char Y = FLOOR(255 * (a - X)); //fractional part from 0 to 255
		switch(X)
		{
		    case 0: red = 0; green = 255-Y; blue = 255 - Y; break;
		    case 1: red = Y; green = Y; blue = Y; break;
		    case 2: red = 255; green = 255-Y; blue = 255 - Y; break;
		    case 3: red = 255; green = 0; blue = 0; break;
		}

		break;
	}
	case RAINBOWSCALE:
	{
		const RFLOAT a = (255 - grey) / 64.; //invert and group
		const int X = FLOOR(a);
		const unsigned char Y = FLOOR(255 * (a - X)); //fractional part from 0 to 255
		switch(X)
		{
		    case 0: red = 255; green = Y; blue = 0; break;
		    case 1: red = 255 - Y; green = 255; blue = 0; break;
		    case 2: red = 0; green = 255; blue = Y; break;
		    case 3: red = 0; green = 255-Y; blue = 255; break;
		    case 4: red = 0; green = 0; blue = 255; break;
		}

		break;
	}
	case CYANBLACKYELLOWSCALE:
	{
		const RFLOAT d_rb = 3 * (grey - 128);
		const RFLOAT d_g = 3 * (std::abs(grey - 128) - 42);
		red   = (unsigned char)(FLOOR(XMIPP_MIN(255., XMIPP_MAX(0.0,  d_rb))));
		green = (unsigned char)(FLOOR(XMIPP_MIN(255., XMIPP_MAX(0.0,  d_g))));
		blue  = (unsigned char)(FLOOR(XMIPP_MIN(255., XMIPP_MAX(0.0, -d_rb))));

		break;
	}
	}
	return;
}
#endif
