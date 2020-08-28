#if defined WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#ifndef _PD_RADAR_H
#define _PD_RADAR_H

#include <fftw3.h>

EXPORT void call_LFM(fftwf_complex *chirp,fftwf_complex *coeff,float timeWidth,float bandWidth,
    float f0);

EXPORT void generate_object(fftwf_complex *chirp,fftwf_complex *signalAll,int targetNumber,
    int *signalPower,int *targetDistance,int *targetVelocity,float RF,int noisePower);

EXPORT void call_PC(fftwf_complex *signalAll,fftwf_complex *coeff,fftwf_complex *pc,bool flag);

EXPORT void call_MTI(fftwf_complex *pc,fftwf_complex *mti);

EXPORT void call_MTD(fftwf_complex *mti,fftwf_complex *mtd,float *freqWeight);

EXPORT void call_CFAR(fftwf_complex *mtd,float *cfar,float *threshold,int protectNumber,
    int refWindowSize,double alpha,double pf,int cfarType,int detectorType);

#endif