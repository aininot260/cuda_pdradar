#ifndef _GPU_H
#define _GPU_H

#include <cstdio>
#include <cufft.h>
#include <cublas_v2.h>

#define CHECK(call)\
{\
  const cudaError_t error=call;\
  if(error!=cudaSuccess)\
  {\
      printf("ERROR: %s:%d,",__FILE__,__LINE__);\
      printf("code:%d,reason:%s\n",error,cudaGetErrorString(error));\
      exit(1);\
  }\
}

extern cufftComplex *pingBuf,*pongBuf;
extern double *window1,*window2;
extern cufftHandle *plan1,*plan2;
extern float *freqWeight_gpu,*cfar_gpu,*threshold_gpu;
extern cublasHandle_t *handle;
extern int streamNum,currentStream;
extern cudaStream_t *stream;

extern float* d_block_sums;
extern float* d_dummy_blocks_sums;
extern float* d_in_block_sums;

float *get_d_block_sums();
float *get_d_in_block_sums();
float *get_d_dummy_blocks_sums();

cufftComplex *get_pong_buf();
cufftComplex *get_ping_buf();
float *get_cfar_gpu();
float *get_threshold_gpu();

#endif