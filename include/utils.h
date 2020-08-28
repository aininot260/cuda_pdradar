#if defined WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#ifndef _UTILS_H
#define _UTILS_H

#ifdef _WIN32
#	include <windows.h>
#else
#	include <sys/time.h>
#endif
#include <cmath>
#include <ctime>
#include <fftw3.h>

#define __C 3.0e8
#define PI 3.14159265358979

#define CFAR_CA   1
#define CFAR_SOCA 2
#define CFAR_GOCA 3

#define Detector_Evelop 4
#define Detector_Square 5
#define Detector_LogSquare 6
#define Detector_Log 7

#define WINDOW_TYPE_RECTANGLE	  1
#define WINDOW_TYPE_BARTLETT	  2
#define WINDOW_TYPE_HANNING		  3
#define WINDOW_TYPE_HAMMING		  4
#define WINDOW_TYPE_BLACKMAN	  5
#define WINDOW_TYPE_STEEPBLACKMAN 6
#define WINDOW_TYPE_KAISER		  7

#ifdef _WIN32
int gettimeofday(struct timeval *tp, void *tzp);
#endif

extern int IO_mode;
extern int blindNumber;
extern int sampleNumber;
extern int totalNumber;
extern int fftSize;
extern int pulseNumber;
extern double parameters1;
extern double parameters2;
extern int windowType1;
extern int windowType2;
extern float sampleRate;
extern double totalTime;
extern double start,_end;
extern int setStreamNum;
extern bool mode0_flag;
extern bool gpu_accelerate;

EXPORT void set_gpu_acceletate(bool flag);

EXPORT void clean_timer();

EXPORT void sync_time();

EXPORT void initialize(double parameters1,double parameters2,int windowType1,int windowType2,
    float *freqWeight,float sampleRate,float wavegate);

EXPORT void release();

EXPORT void set_IO_mode(int flag);

EXPORT void set_stream_number(int num);

EXPORT int set_blindNumber(float sampleRate,float timeWidth);

EXPORT int set_sampleNumber(float sampleRate,float PRT);

EXPORT int set_totalNumber(int sampleNumber,int pulseNumber);

EXPORT void get_total_time();

EXPORT void malloc_IO_memory(fftwf_complex **signalAll,fftwf_complex **coeff,float **cfar,float **threshold);

EXPORT void free_IO_memory(fftwf_complex *signalAll,fftwf_complex *coeff,float *cfar,float *threshold);

double cpuSecond();

void create_window(double*, int, int, double);

inline int nextpow2(int n)
{
	double a,b;
	a=log(double(n))/log(2.0);
	b=double(int(a));
	if(fabs(a-b)<1e-9)
		return n;
	else
		return int(pow(2, b + 1));
}

inline double alpha_factorial(int n, int m)
{
	double f=1;
	for(int i=n;i>m;i--) f=f*i;
	return f;
}

#endif