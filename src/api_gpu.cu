#include "utils.h"
#include "utils_gpu.h"
#include "api_gpu.h"

inline int getStreamNum(size_t fftSize)
{
    size_t avail;
    size_t total;
    cudaMemGetInfo(&avail,&total);
    size_t ret=(avail-16*fftSize)/(40*fftSize+10*1024*1024);
    return (int)ret;
}

void init_gpu(int fftSize,int pulseNumber,int sampleNumber,double parameters1,double parameters2,
	int windowType1,int windowType2,float *freqWeight)
{
    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp,dev));
    printf("Using device %d: %s\n",dev,deviceProp.name);
    CHECK(cudaSetDevice(dev));

    streamNum=getStreamNum(fftSize);
    if(streamNum==0)
    {
        puts("ERROR:The GPU memory is not enough.");
        exit(0);
    }
    if(setStreamNum!=0)
    {
        if(setStreamNum>streamNum)
        {
            puts("ERROR:The Manual Stream Number is Too Much.");
            exit(0);
        }
        streamNum=setStreamNum;
    }
    printf("GPU streamNum : %d\n",streamNum);
    stream=(cudaStream_t *)fftw_malloc(streamNum*sizeof(cudaStream_t));
    for(int i=0;i<streamNum;i++)
        cudaStreamCreateWithFlags(&stream[i],cudaStreamDefault/*cudaStreamNonBlocking*/);

    CHECK(cudaMalloc((void**)&pingBuf,sizeof(cufftComplex)*fftSize*streamNum));
    CHECK(cudaMemset(pingBuf,0,sizeof(cufftComplex)*fftSize*streamNum));
    CHECK(cudaMalloc((void**)&pongBuf,sizeof(cufftComplex)*fftSize*streamNum));
    CHECK(cudaMemset(pongBuf,0,sizeof(cufftComplex)*fftSize*streamNum));
    CHECK(cudaMalloc((void**)&cfar_gpu,sizeof(float)*pulseNumber*sampleNumber*streamNum));
    CHECK(cudaMalloc((void**)&threshold_gpu,sizeof(float)*pulseNumber*sampleNumber*streamNum));

    CHECK(cudaMalloc((void**)&d_block_sums,sizeof(float)*(fftSize/1024+1)*streamNum));
    // CHECK(cudaMemset(d_block_sums,0,sizeof(float)*(fftSize/1024+1)*streamNum));
    CHECK(cudaMalloc((void**)&d_in_block_sums,sizeof(float)*(fftSize/1024+1)*streamNum));
    // CHECK(cudaMemset(d_in_block_sums,0,sizeof(float)*(fftSize/1024+1)*streamNum));
    CHECK(cudaMalloc((void**)&d_dummy_blocks_sums,sizeof(float)*streamNum));
    // CHECK(cudaMemset(d_dummy_blocks_sums,0,sizeof(float)*streamNum));


    plan1=(cufftHandle *)fftw_malloc(sizeof(cufftHandle)*streamNum);
    plan2=(cufftHandle *)fftw_malloc(sizeof(cufftHandle)*streamNum);
    handle=(cublasHandle_t *)fftw_malloc(sizeof(cublasHandle_t)*streamNum);
    for(int i=0;i<streamNum;i++)
    {
        cufftPlan1d(&plan1[i],fftSize,CUFFT_C2C,1);
        cufftPlan1d(&plan2[i],pulseNumber-1,CUFFT_C2C,sampleNumber);
        cublasCreate(&handle[i]);
        cufftSetStream(plan1[i],stream[i]);
        cufftSetStream(plan2[i],stream[i]);
        cublasSetStream(handle[i],stream[i]);
    }

    CHECK(cudaMalloc((void**)&window1,sizeof(double)*fftSize));
    CHECK(cudaMalloc((void**)&window2,sizeof(double)*pulseNumber));
    double *window_func=(double *)fftw_malloc(fftSize*sizeof(double));
    create_window(window_func,fftSize,windowType1,parameters1);
    CHECK(cudaMemcpy(window1,window_func,sizeof(double)*fftSize,cudaMemcpyHostToDevice));
    fftw_free(window_func);
    window_func=(double *)fftw_malloc(pulseNumber*sizeof(double));
    create_window(window_func,pulseNumber,windowType2,parameters2);
    CHECK(cudaMemcpy(window2,window_func,sizeof(double)*pulseNumber,cudaMemcpyHostToDevice));
    fftw_free(window_func);
    CHECK(cudaMalloc((void**)&freqWeight_gpu,sizeof(float)*pulseNumber*sampleNumber));
    CHECK(cudaMemcpy(freqWeight_gpu,freqWeight,sizeof(float)*pulseNumber*sampleNumber,cudaMemcpyHostToDevice));
}

void release_gpu()
{
    CHECK(cudaFree(pingBuf));
    CHECK(cudaFree(pongBuf));
    CHECK(cudaFree(cfar_gpu));
    CHECK(cudaFree(threshold_gpu));

    CHECK(cudaFree(d_block_sums));
    CHECK(cudaFree(d_in_block_sums));
    CHECK(cudaFree(d_dummy_blocks_sums));

    for(int i=0;i<streamNum;i++)
    {
        cudaStreamDestroy(stream[i]);
        cufftDestroy(plan1[i]);
        cufftDestroy(plan2[i]);
        cublasDestroy(handle[i]);
    }

    CHECK(cudaFree(window1));
    CHECK(cudaFree(window2));
    CHECK(cudaFree(freqWeight_gpu));

    fftw_free(stream);
    fftw_free(plan1);
    fftw_free(plan2);
    fftw_free(handle);
}

void malloc_IO_memory_gpu(fftwf_complex **signalAll,fftwf_complex **coeff,float **cfar,float **threshold)
{
    cudaMallocHost((fftwf_complex**)signalAll,totalNumber*sizeof(fftwf_complex));
    cudaMallocHost((fftwf_complex**)coeff,blindNumber*sizeof(fftwf_complex));
    cudaMallocHost((float**)cfar,(totalNumber)*sizeof(float));
    cudaMallocHost((float**)threshold,(totalNumber)*sizeof(float));
}

void free_IO_memory_gpu(fftwf_complex *signalAll,fftwf_complex *coeff,float *cfar,float *threshold)
{
    cudaFreeHost(coeff);
    cudaFreeHost(signalAll);
    cudaFreeHost(cfar);
    cudaFreeHost(threshold);
}

void sync_time_gpu()
{
    cudaDeviceSynchronize();
    _end=cpuSecond();
    double timeCost=_end-start;
	totalTime+=timeCost;
}