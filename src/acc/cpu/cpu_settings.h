#ifndef CPU_SETTINGS_H_
#define CPU_SETTINGS_H_

#include "src/acc/settings.h"


// GENERAL -----------------------------
#define MAX_RESOL_SHARED_MEM 		32
#define BLOCK_SIZE  				128
// -------------------------------------


// COARSE DIFF -------------------------
#define PREFETCH_FRACTION_3D 4  //4
#define PREFETCH_FRACTION_2D 2  //2

#define D2C_BLOCK_SIZE_2D           256	//256
#define D2C_EULERS_PER_BLOCK_2D     16	//16
#define D2C_BLOCK_SIZE_REF3D        256 //256
#define D2C_EULERS_PER_BLOCK_REF3D  16	//16
#define D2C_BLOCK_SIZE_DATA3D       64	//64
#define D2C_EULERS_PER_BLOCK_DATA3D 32	//32

#define D2C_CC_BLOCK_SIZE_2D           256	//256
#define D2C_CC_EULERS_PER_BLOCK_2D     16	//16
#define D2C_CC_BLOCK_SIZE_REF3D        256 //256
#define D2C_CC_EULERS_PER_BLOCK_REF3D  16	//16
#define D2C_CC_BLOCK_SIZE_DATA3D       64	//64
#define D2C_CC_EULERS_PER_BLOCK_DATA3D 32	//32
// -------------------------------------


// FINE DIFF ---------------------------
#define D2F_BLOCK_SIZE_2D 			8	//8
#define D2F_CC_BLOCK_SIZE_2D 		8	//8
#define D2F_CHUNK_2D 				7	//7

#define D2F_BLOCK_SIZE_REF3D 		8	//8
#define D2F_CC_BLOCK_SIZE_REF3D 	8	//8
#define D2F_CHUNK_REF3D 			7	//7

#define D2F_BLOCK_SIZE_DATA3D 		8	//8
#define D2F_CC_BLOCK_SIZE_DATA3D 	8	//8
#define D2F_CHUNK_DATA3D			4	//4
// -------------------------------------


// WAVG --------------------------------
#define WAVG_BLOCK_SIZE_DATA3D	 	512	//512
#define WAVG_BLOCK_SIZE_REF3D	 	512	//512
#define WAVG_BLOCK_SIZE_2D 	 		256	//256
// -------------------------------------


// MISC --------------------------------
#define SUMW_BLOCK_SIZE 	  	32
#define SOFTMASK_BLOCK_SIZE 	128
#define CFTT_BLOCK_SIZE 	 	128
#define PROBRATIO_BLOCK_SIZE 	128
#define POWERCLASS_BLOCK_SIZE 	128
#define PROJDIFF_CHUNK_SIZE 	14

// -------------------------------------

// RANDOMIZATION -----------------------
#define RND_BLOCK_NUM                   64
#define RND_BLOCK_SIZE                  32
// -------------------------------------


#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_PREFETCH_COUNT 3
#define BP_2D_BLOCK_SIZE		128	//128
#define BP_REF3D_BLOCK_SIZE		128	//128
#define BP_DATA3D_BLOCK_SIZE	512	//640

#define REF_GROUP_SIZE 3			// -- Number of references to be treated per block --
									// This applies to wavg and reduces global memory
									// accesses roughly proportionally, but scales shared
									// memory usage by allocating
									// ( 6*REF_GROUP_SIZE + 4 ) * BLOCK_SIZE XFLOATS. // DEPRECATED

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
 #ifndef USE_INTEL_COMPILER
  #define USE_INTEL_COMPILER
 #endif
 #include <mathimf.h>
#endif
#define RESTRICT __restrict

// This is for AVX512
#ifdef ACC_DOUBLE_PRECISION
 #define SIMD_LEN 8
#else
 #define SIMD_LEN 16
#endif

#ifndef MEM_ALIGN
 #define MEM_ALIGN 64
#endif


#endif /* CPU_SETTINGS_H_ */
