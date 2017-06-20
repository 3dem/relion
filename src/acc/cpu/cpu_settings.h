#ifndef CPU_SETTINGS_H_
#define CPU_SETTINGS_H_

#include "src/acc/settings.h"


// GENERAL -----------------------------
#define MAX_RESOL_SHARED_MEM 		32
#define BLOCK_SIZE  				128
// -------------------------------------


// COARSE DIFF -------------------------
#define D2C_BLOCK_SIZE_2D 			512
#define D2C_EULERS_PER_BLOCK_2D 	4

#define D2C_BLOCK_SIZE_REF3D 		128
#define D2C_EULERS_PER_BLOCK_REF3D 	16

#define D2C_BLOCK_SIZE_DATA3D 		64
#define D2C_EULERS_PER_BLOCK_DATA3D 32
// -------------------------------------


// FINE DIFF ---------------------------
#define D2F_BLOCK_SIZE_2D 			8
#define D2F_CHUNK_2D 				7

#define D2F_BLOCK_SIZE_REF3D 		8
#define D2F_CHUNK_REF3D 			7

#define D2F_BLOCK_SIZE_DATA3D 		8
#define D2F_CHUNK_DATA3D			4
// -------------------------------------


// WAVG --------------------------------
#define WAVG_BLOCK_SIZE_DATA3D	 	512
#define WAVG_BLOCK_SIZE 	 		256
// -------------------------------------


// MISC --------------------------------
#define SUMW_BLOCK_SIZE 	  	32
#define SOFTMASK_BLOCK_SIZE 	128
#define CFTT_BLOCK_SIZE 	 	128
#define PROBRATIO_BLOCK_SIZE 	128
#define POWERCLASS_BLOCK_SIZE 	128
#define PROJDIFF_CHUNK_SIZE 	14

#define MEM_ALIGN               64

#ifdef __INTEL_COMPILER
#  include <mathimf.h>
#  define RESTRICT __restrict__
#else
#  define RESTRICT
#endif


#endif /* CPU_SETTINGS_H_ */
