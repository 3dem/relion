#ifndef ACC_SETTINGS_H_
#define ACC_SETTINGS_H_

#ifdef ACC_DOUBLE_PRECISION
	#define XFLOAT double
	typedef struct{ XFLOAT x; XFLOAT y;} double2;
	#define ACCCOMPLEX double2
#else
	#define XFLOAT float
    typedef struct{ XFLOAT x; XFLOAT y;} float2;
	#define ACCCOMPLEX float2
#endif

#ifdef RELION_SINGLE_PRECISION
	#define RFLOAT float
#else
	#define RFLOAT double
#endif

#endif /* ACC_SETTINGS_H_ */
