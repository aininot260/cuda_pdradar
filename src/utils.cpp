#include <cstdio>

#ifdef USE_GPU
#include "api_gpu.h"
#endif
#include "utils.h"

int IO_mode=1;
int blindNumber;
int sampleNumber;
int totalNumber;
int fftSize;
int pulseNumber;
double parameters1;
double parameters2;
int windowType1;
int windowType2;
float sampleRate;
double totalTime;
double start=0.0,_end=0.0;
int setStreamNum=0;
bool mode0_flag=false;
bool gpu_accelerate=true;

void set_gpu_acceletate(bool flag)
{
    gpu_accelerate=flag;
}

void clean_timer()
{
    mode0_flag=0;
    totalTime=0.0;
}

void sync_time()
{
    if(IO_mode) return;
    if(gpu_accelerate)
    {
        #ifdef USE_GPU
        sync_time_gpu();
        #endif
    }
}

void initialize(double _parameters1,double _parameters2,int _windowType1,int _windowType2,
    float *freqWeight,float _sampleRate,float wavegate)
{
    double start,_end;
    start=cpuSecond();

    parameters1=_parameters1;
    parameters2=_parameters2;
    windowType1=_windowType1;
    windowType2=_windowType2;
    sampleRate=_sampleRate;
    fftSize=nextpow2(sampleRate*wavegate)*pulseNumber;
    printf("fftSize : %d\n",fftSize);

    #ifdef USE_GPU
    init_gpu(fftSize,pulseNumber,sampleNumber,parameters1,parameters2,windowType1,windowType2,freqWeight);
    #endif

    _end=cpuSecond();
    double timeCost=_end-start;
	printf("Initialize Time Cost: %.6lfs\n",timeCost);
}

void release()
{
    double start,_end;
    start=cpuSecond();

    #ifdef USE_GPU
    release_gpu();
    #endif
    
    _end=cpuSecond();
    double timeCost=_end-start;
	printf("Release Time Cost: %.6lfs\n",timeCost);
}

void set_IO_mode(int flag)
{
    IO_mode=flag;
    printf("IO mode : %d\n",IO_mode);
}

void set_stream_number(int num)
{
    #ifndef USE_GPU
    puts("WARNING:This Version of PD_Radar Don't Support GPU ACCELERATED.");
    return;
    #endif
    setStreamNum=num;
    puts("The Manual Stream mode is ON.");
}

int set_blindNumber(float sampleRate,float timeWidth)
{
    blindNumber=floor(sampleRate*timeWidth); //回波的采样点数
    if(blindNumber%2!=0) blindNumber++;
    printf("blindNumber : %d\n",blindNumber);
    return blindNumber;
}

int set_sampleNumber(float sampleRate,float PRT)
{
    sampleNumber=floor(sampleRate*PRT);
    printf("sampleNumber : %d\n",sampleNumber);
    return sampleNumber;
}

int set_totalNumber(int sampleNumber,int _pulseNumber)
{
    totalNumber=sampleNumber*_pulseNumber;
    pulseNumber=_pulseNumber;
    printf("pulseNumber : %d\n",pulseNumber);
    printf("totalNumber : %d\n",totalNumber);
    return totalNumber;
}

double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return((double)tp.tv_sec+(double)tp.tv_usec*1e-6);
}

void get_total_time()
{
    if(IO_mode) return;
    printf("Total Signal Produce Time Cost: %.6lfs\n",totalTime);
}

void malloc_IO_memory(fftwf_complex **signalAll,fftwf_complex **coeff,float **cfar,float **threshold)
{
    #ifdef USE_GPU
    malloc_IO_memory_gpu(signalAll,coeff,cfar,threshold);
    #else
    *coeff=(fftwf_complex *)fftw_malloc(blindNumber*sizeof(fftwf_complex));
    *signalAll=(fftwf_complex *)fftw_malloc(totalNumber*sizeof(fftwf_complex));
    *cfar=(float *)fftw_malloc((totalNumber)*sizeof(float));
    *threshold=(float *)fftw_malloc((totalNumber)*sizeof(float));
    #endif
}

void free_IO_memory(fftwf_complex *signalAll,fftwf_complex *coeff,float *cfar,float *threshold)
{
    #ifdef USE_GPU
    free_IO_memory_gpu(signalAll,coeff,cfar,threshold);
    #else
    fftw_free(coeff);
    fftw_free(signalAll);
    fftw_free(cfar);
    fftw_free(threshold);
    #endif
}

#ifdef _WIN32
int gettimeofday(struct timeval *tp, void *tzp)
{
    time_t clock;
    struct tm tm;
    SYSTEMTIME wtm;
    GetLocalTime(&wtm);
    tm.tm_year   = wtm.wYear - 1900;
    tm.tm_mon   = wtm.wMonth - 1;
    tm.tm_mday   = wtm.wDay;
    tm.tm_hour   = wtm.wHour;
    tm.tm_min   = wtm.wMinute;
    tm.tm_sec   = wtm.wSecond;
    tm. tm_isdst  = -1;
    clock = mktime(&tm);
    tp->tv_sec = clock;
    tp->tv_usec = wtm.wMilliseconds * 1000;
    return (0);
}
#endif