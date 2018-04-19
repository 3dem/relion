/******************************************************************************
 **        Title: matlab.h
 **  Description: Read and write tArrays to/from matlab files.
 **
 **       Author: Brian Amberg, 2007
 **               Computer Science Department, University Basel (CH)
 **
 ** Linking:
 **  When using these functions you should link with the following commandline:
 **
 **  -I$(MATLAB)/extern/include -L$(MATLAB)/bin/em64 -L$(MATLAB)/bin/glnx86 \
 **    -lmat -lmx -lut -licui18n -licuio -licuuc -licudata -lhdf5 \
 **    -D__GRAVIS__MATLAB__
 **
 **  And make shure that $(MATLAB)/bin/glnx86 is in your LD_LIBRARY_PATH
 **
 ******************************************************************************/
#ifndef __GRAVIS__IO__ARRAY__MATLAB__
#define __GRAVIS__IO__ARRAY__MATLAB__

#include <gravis/tArray.h>
#include <gravis/Exception.h>

#ifdef __GRAVIS__MATLAB__
#include <mat.h>
#include <matrix.h>
#include <string>
#include "../../StringFormat.h"
#include "../../tRGB.h"
#include "../../tRGBA.h"
#include "../../t2Vector.h"
#include "../../t3Vector.h"
#include "../../t4Vector.h"
#include "../../Tuple.h"
#include "../matlabFileIO.h"

namespace gravis
{
  namespace io
  {
    class ArrayMatlab: public MatlabFileIO
    {
      private:
        // Setup knowledge about how to make an tArray out of this mxArray type
#define DEFINE_MX_ARRAY_CONVERSION(type, classid) \
	static inline void fillArray(tArray<type>               &a, const mxArray *mx)                                                                                                          \
	  { const size_t channels=sizeof(a[0])/sizeof(type);                                                                                                                      \
	    if (mxGetClassID(mx) != classid) GRAVIS_THROW3(gravis::Exception, "Incompatible datatype", StringFormat("Found ")(mxGetClassName(mx))(" and expected ")( #classid )); \
	    if (mxGetM(mx) != channels) GRAVIS_THROW3(gravis::Exception, "Array has the wrong number of rows.", StringFormat(mxGetM(mx))(" found and ")(channels)(" expected"));  \
	    a.resize(mxGetN(mx)); memcpy(&a[0], mxGetPr(mx), channels*sizeof(type) * a.size());    }                                                                              \
	static inline void fillArray(tArray< t2Vector< type > > &a, const mxArray *mx)                                                                                            \
	  { const size_t channels=sizeof(a[0])/sizeof(type);                                                                                                                      \
	    if (mxGetClassID(mx) != classid) GRAVIS_THROW3(gravis::Exception, "Incompatible datatype", StringFormat("Found ")(mxGetClassName(mx))(" and expected ")( #classid )); \
	    if (mxGetM(mx) != channels) GRAVIS_THROW3(gravis::Exception, "Array has the wrong number of rows.", StringFormat(mxGetM(mx))(" found and ")(channels)(" expected"));  \
	    a.resize(mxGetN(mx)); memcpy(&a[0], mxGetPr(mx), channels*sizeof(type) * a.size());    }                                                                              \
	static inline void fillArray(tArray< t3Vector< type > > &a, const mxArray *mx)                                                                                            \
	  { const size_t channels=sizeof(a[0])/sizeof(type);                                                                                                                      \
	    if (mxGetClassID(mx) != classid) GRAVIS_THROW3(gravis::Exception, "Incompatible datatype", StringFormat("Found ")(mxGetClassName(mx))(" and expected ")( #classid )); \
	    if (mxGetM(mx) != channels) GRAVIS_THROW3(gravis::Exception, "Array has the wrong number of rows.", StringFormat(mxGetM(mx))(" found and ")(channels)(" expected"));  \
	    a.resize(mxGetN(mx)); memcpy(&a[0], mxGetPr(mx), channels*sizeof(type) * a.size());    }                                                                              \
	static inline void fillArray(tArray< t4Vector< type > > &a, const mxArray *mx)                                                                                            \
	  { const size_t channels=sizeof(a[0])/sizeof(type);                                                                                                                      \
	    if (mxGetClassID(mx) != classid) GRAVIS_THROW3(gravis::Exception, "Incompatible datatype", StringFormat("Found ")(mxGetClassName(mx))(" and expected ")( #classid )); \
	    if (mxGetM(mx) != channels) GRAVIS_THROW3(gravis::Exception, "Array has the wrong number of rows.", StringFormat(mxGetM(mx))(" found and ")(channels)(" expected"));  \
	    a.resize(mxGetN(mx)); memcpy(&a[0], mxGetPr(mx), channels*sizeof(type) * a.size());    }                                                                              \
	static inline void fillArray(tArray< tRGB< type > >     &a, const mxArray *mx)                                                                                            \
	  { const size_t channels=sizeof(a[0])/sizeof(type);                                                                                                                      \
	    if (mxGetClassID(mx) != classid) GRAVIS_THROW3(gravis::Exception, "Incompatible datatype", StringFormat("Found ")(mxGetClassName(mx))(" and expected ")( #classid )); \
	    if (mxGetM(mx) != channels) GRAVIS_THROW3(gravis::Exception, "Array has the wrong number of rows.", StringFormat(mxGetM(mx))(" found and ")(channels)(" expected"));  \
	    a.resize(mxGetN(mx)); memcpy(&a[0], mxGetPr(mx), channels*sizeof(type) * a.size());    }                                                                              \
	static inline void fillArray(tArray< tRGBA< type > >    &a, const mxArray *mx)                                                                                            \
	  { const size_t channels=sizeof(a[0])/sizeof(type);                                                                                                                      \
	    if (mxGetClassID(mx) != classid) GRAVIS_THROW3(gravis::Exception, "Incompatible datatype", StringFormat("Found ")(mxGetClassName(mx))(" and expected ")( #classid )); \
	    if (mxGetM(mx) != channels) GRAVIS_THROW3(gravis::Exception, "Array has the wrong number of rows.", StringFormat(mxGetM(mx))(" found and ")(channels)(" expected"));  \
	    a.resize(mxGetN(mx)); memcpy(&a[0], mxGetPr(mx), channels*sizeof(type) * a.size());    }

        DEFINE_MX_ARRAY_CONVERSION(float,   mxSINGLE_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(double,   mxDOUBLE_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int8_t,   mxINT8_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int16_t,  mxINT16_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int32_t,  mxINT32_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int64_t,  mxINT64_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint8_t,  mxUINT8_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint16_t, mxUINT16_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint32_t, mxUINT32_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint64_t, mxUINT64_CLASS);
#undef DEFINE_MX_ARRAY_CONVERSION

        // Setup knowledge about how to make an mxArray out of the gravis datatypes
#define DEFINE_MX_ARRAY_CONVERSION(type, classid) \
	static inline mxArray *mxFromArray(const tArray<type> &a)                 { mxArray *mx = mxCreateNumericMatrix(1, a.size(), classid, mxREAL); memcpy(mxGetPr(mx), &a[0], sizeof(type) * a.size());   return mx; } \
	static inline mxArray *mxFromArray(const tArray< t2Vector< type > > &a)   { mxArray *mx = mxCreateNumericMatrix(2, a.size(), classid, mxREAL); memcpy(mxGetPr(mx), &a[0], 2*sizeof(type) * a.size()); return mx; } \
	static inline mxArray *mxFromArray(const tArray< t3Vector< type > > &a)   { mxArray *mx = mxCreateNumericMatrix(3, a.size(), classid, mxREAL); memcpy(mxGetPr(mx), &a[0], 3*sizeof(type) * a.size()); return mx; } \
	static inline mxArray *mxFromArray(const tArray< t4Vector< type > > &a)   { mxArray *mx = mxCreateNumericMatrix(4, a.size(), classid, mxREAL); memcpy(mxGetPr(mx), &a[0], 4*sizeof(type) * a.size()); return mx; } \
	static inline mxArray *mxFromArray(const tArray< tRGB< type > > &a)       { mxArray *mx = mxCreateNumericMatrix(3, a.size(), classid, mxREAL); memcpy(mxGetPr(mx), &a[0], 3*sizeof(type) * a.size()); return mx; } \
	static inline mxArray *mxFromArray(const tArray< tRGBA< type > > &a)      { mxArray *mx = mxCreateNumericMatrix(4, a.size(), classid, mxREAL); memcpy(mxGetPr(mx), &a[0], 4*sizeof(type) * a.size()); return mx; }

        static inline mxArray* mxFromArray(const tArray< Tuple2 > &a)
        {
          mxArray* mx = mxCreateNumericMatrix(2, a.size(), mxINT32_CLASS, mxREAL);
          memcpy(mxGetPr(mx), &a[0], 2*sizeof(int) * a.size());
          return mx;
        }
        static inline mxArray* mxFromArray(const tArray< Tuple3 > &a)
        {
          mxArray* mx = mxCreateNumericMatrix(3, a.size(), mxINT32_CLASS, mxREAL);
          memcpy(mxGetPr(mx), &a[0], 3*sizeof(int) * a.size());
          return mx;
        }

        DEFINE_MX_ARRAY_CONVERSION(float,   mxSINGLE_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(double,  mxDOUBLE_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int8_t,  mxINT8_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int16_t, mxINT16_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int32_t, mxINT32_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(int64_t, mxINT64_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint8_t,  mxUINT8_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint16_t, mxUINT16_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint32_t, mxUINT32_CLASS);
        DEFINE_MX_ARRAY_CONVERSION(uint64_t, mxUINT64_CLASS);
#undef DEFINE_MX_ARRAY_CONVERSION

      public:
        using MatlabFileIO::get;
        using MatlabFileIO::put;
        using MatlabFileIO::hasVar;
        using MatlabFileIO::getNDim;
        using MatlabFileIO::getDimensions;
        using MatlabFileIO::getClassID;

        ArrayMatlab():MatlabFileIO() {}
        ArrayMatlab(const std::string& filename, const std::string& mode):MatlabFileIO(filename, mode) {};

        /**
         * Load an array from an matlab file.
         **/
        template<class T>
        void get(tArray<T> &out, const std::string& varname = "gravis_array")
        {
          if(!pmat)
            GRAVIS_THROW2(gravis::Exception, "Error file not open! Call open first.");

          /*
           * Read in each array we just wrote
           */
          mxArray* pa = matGetVariable(pmat, varname.c_str());
          if (pa == NULL)
          {
            matClose(pmat);
            GRAVIS_THROW3(gravis::Exception, "Did not find variable in file.", varname);
          }

          fillArray(out, pa);

          /* clean up before exit */
          mxDestroyArray(pa);
        }

        /**
         * Load an array from an matlab file.
         **/
        template<class T>
        static void read(tArray<T> &out, const std::string& filename, const std::string& varname = "gravis_array")
        {
          ArrayMatlab am(filename,"r");
          am.get(out,varname);
        }


        /**
         * Save an array to a matlab file. This should work with float and double data.	 **/
        template<class T>
        void put(const tArray<T> &in, const std::string& varname = "gravis_array")
        {
          if(!pmat)
            GRAVIS_THROW2(gravis::Exception, "Error file not open! Call open first.");

          mxArray* pa = mxFromArray(in);
          if (pa == NULL) GRAVIS_THROW3(gravis::Exception, "Could not convert to mxArray. Can not save: ", varname);
          int status;

          status = matPutVariable(pmat, varname.c_str(), pa);
          if (status != 0)
          {
            mxDestroyArray(pa);
            matClose(pmat);
            GRAVIS_THROW3(gravis::Exception, "Could not put variable into file. Matlab error code: ", StringFormat(status));
          }
          mxDestroyArray(pa);
        }

        /**
         * Save an array to a matlab file. This should work with float and double data.	 **/
        template<class T>
        static void write(const std::string& filename, const tArray<T> &in, const std::string& varname = "gravis_array")
        {
          ArrayMatlab am(filename,"w");
          am.put(in,varname);
        }
    };
  }

} // end namespace

#endif
#endif
