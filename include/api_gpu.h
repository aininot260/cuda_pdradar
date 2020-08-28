#ifndef _API_GPU_H
#define _API_GPU_H

#include <fftw3.h>

void init_gpu(int fftSize,int pulseNumber,int sampleNumber,double parameters1,double parameters2,
	int windowType1,int windowType2,float *freqWeight);

void release_gpu();

void call_pc_gpu(fftwf_complex *signalAll,fftwf_complex *coeff,fftwf_complex *pc,int totalNumber,int blindNumber,
    int fftSize,bool flag,double parameters,int windowType);

void call_mti_gpu(fftwf_complex *pc,fftwf_complex *mti,int &pulseNumber,int sampleNumber,int &totalNumber);

void call_mtd_gpu(fftwf_complex *mti,fftwf_complex *mtd,int sampleNumber,int pulseNumber,int totalNumber,
    int windowType,double parameters,float *freqWeight);

void call_cfar_gpu(fftwf_complex *mtd,float *cfar,float *threshold,int totalNumber,int pulseNumber,int sampleNumber,
    int protectNumber,int refWindowSize,double alpha,double pf,int cfarType,int detectorType,double reference);

void malloc_IO_memory_gpu(fftwf_complex **signalAll,fftwf_complex **coeff,float **cfar,float **threshold);

void free_IO_memory_gpu(fftwf_complex *signalAll,fftwf_complex *coeff,float *cfar,float *threshold);

void sync_time_gpu();

#endif