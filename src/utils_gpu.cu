#include "utils.h"
#include "utils_gpu.h"

cufftComplex *pingBuf,*pongBuf;
double *window1,*window2;
float *freqWeight_gpu,*cfar_gpu,*threshold_gpu;
cufftHandle *plan1,*plan2;
cublasHandle_t *handle;
int streamNum,currentStream=-1;
cudaStream_t *stream;

float* d_block_sums;
float* d_dummy_blocks_sums;
float* d_in_block_sums;

float *get_d_block_sums()
{
    return d_block_sums+(fftSize/1024+1)*currentStream;
}

float *get_d_in_block_sums()
{
    return d_in_block_sums+(fftSize/1024+1)*currentStream;
}

float *get_d_dummy_blocks_sums()
{
    return d_dummy_blocks_sums+1*currentStream;
}

cufftComplex *get_pong_buf()
{
    return pongBuf+currentStream*fftSize;
}

cufftComplex *get_ping_buf()
{
    return pingBuf+currentStream*fftSize;
}

float *get_cfar_gpu()
{
    return cfar_gpu+pulseNumber*sampleNumber*currentStream;
}
float *get_threshold_gpu()
{
    return threshold_gpu+pulseNumber*sampleNumber*currentStream;
}