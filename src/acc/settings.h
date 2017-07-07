#ifndef ACC_SETTINGS_H_
#define ACC_SETTINGS_H_

#ifdef ACC_DOUBLE_PRECISION
	#define XFLOAT double
	#ifndef CUDA
		typedef struct{ XFLOAT x; XFLOAT y;} double2;
	#endif
	#define ACCCOMPLEX double2
#else
	#define XFLOAT float
	#ifndef CUDA
		typedef struct{ XFLOAT x; XFLOAT y;} float2;
	#endif
	#define ACCCOMPLEX float2
#endif

#ifdef RELION_SINGLE_PRECISION
	#define RFLOAT float
#else
	#define RFLOAT double
#endif

#endif /* ACC_SETTINGS_H_ */
