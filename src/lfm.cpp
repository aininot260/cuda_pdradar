#include <cstdio>
#include <cmath>
#include <cstring>
#include <assert.h> 
#include <algorithm>
#include <cstdlib>

#include "pd_radar.h"
#include "utils.h"

using namespace std;

void call_LFM(fftwf_complex *chirp,fftwf_complex *coeff,float timeWidth,float bandWidth,float f0)
{
    start=cpuSecond();

    float mu=bandWidth/timeWidth;
    for(int i=0;i<blindNumber;i++)
    {
        float t=i*1.0/sampleRate-timeWidth/2;
        chirp[i][0]=cos(2*PI*f0*i*1.0/sampleRate+PI*mu*t*t);
        chirp[i][1]=sin(2*PI*f0*i*1.0/sampleRate+PI*mu*t*t);
    }
    for(int i=0;i<blindNumber;i++)
    {
        coeff[i][0]=chirp[blindNumber-i-1][0];
        coeff[i][1]=-chirp[blindNumber-i-1][1];
    }
    
    _end=cpuSecond();
    double timeCost=_end-start;
	printf("Call LFM Time Cost: %.6lfs\n",timeCost);

    if(IO_mode==1)
    {
        FILE *fp=fopen("LFM.dat", "w");
        for (int i=0;i<blindNumber;i++)
            fprintf(fp,"%f\n%f\n",chirp[i][0],chirp[i][1]);
        fclose(fp);
    }
}

inline void normrnd(float avg,float sigma,fftwf_complex *systemNoise,int totalNumber)
{
    srand((unsigned)time(0));

    for(int i=0;i<totalNumber;i++)
    {
        float x,dScope,y,fengzhi;
        do{
            x=((float)rand()/RAND_MAX)*6*sigma+avg-3*sigma;
            y = 1.0/(sqrt(2*PI)*sigma)*exp(-1*(x-avg)*(x-avg)/(2*sigma*sigma));
            fengzhi = 1.0/(sqrt(2*PI)*sigma);
            dScope = ((float)rand()/RAND_MAX)*fengzhi;
        }while(dScope>y);
        systemNoise[i][0]=x;
        do{
            x=((float)rand()/RAND_MAX)*6*sigma+avg-3*sigma;
            y = 1.0/(sqrt(2*PI)*sigma)*exp(-1*(x-avg)*(x-avg)/(2*sigma*sigma));
            fengzhi = 1.0/(sqrt(2*PI)*sigma);
            dScope = ((float)rand()/RAND_MAX)*fengzhi;
        }while(dScope>y);
        systemNoise[i][1]=x;
    }
}

void generate_object(fftwf_complex *chirp,fftwf_complex *signalAll,int targetNumber,int *signalPower,
    int *targetDistance,int *targetVelocity,float RF,int noisePower)
{
    float lambda=__C/RF; //雷达工作波长 电磁波波长
    int *delayNumber=(int *)fftw_malloc(targetNumber*sizeof(int)); //把目标距离换算成采样点（距离门）
    float *targetFd=(float *)fftw_malloc(targetNumber*sizeof(float)); //目标多卜勒
    for(int i=0;i<targetNumber;i++)
        delayNumber[i]=floor(sampleRate*2*targetDistance[i]/__C);

    int maxDelayNumber=0;
    for(int i=0;i<targetNumber;i++)
        maxDelayNumber=max(delayNumber[i],maxDelayNumber);
    
    if(sampleNumber<=blindNumber+maxDelayNumber)
    {
        puts("ERROR:Please insure the sampleNumber > blindNumber + maxDelayNumber.");
        exit(0);
    }

    for(int i=0;i<targetNumber;i++)
        targetFd[i]=2*targetVelocity[i]/lambda;

    fftwf_complex *signalTemp=(fftwf_complex *)fftw_malloc(sampleNumber*sizeof(fftwf_complex)); //一个脉冲
    fftwf_complex *signal=(fftwf_complex *)fftw_malloc(totalNumber*sizeof(fftwf_complex));
    fftwf_complex *freqMove=(fftwf_complex *)fftw_malloc(totalNumber*sizeof(fftwf_complex));
    fftwf_complex *systemNoise=(fftwf_complex *)fftw_malloc(totalNumber*sizeof(fftwf_complex));

    start=cpuSecond();

    for(int i=0;i<targetNumber;i++) //依次产生各个目标
    {
        for(int j=0;j<sampleNumber;j++)
            signalTemp[j][0]=signalTemp[j][1]=0.0;

        for(int j=delayNumber[i];j<delayNumber[i]+blindNumber;j++)
        {
            signalTemp[j][0]=sqrt(signalPower[i])*chirp[j-delayNumber[i]][0];//一个脉冲的1个目标（未加多普勒速度）
            signalTemp[j][1]=sqrt(signalPower[i])*chirp[j-delayNumber[i]][1];
        }

        for(int j=0;j<pulseNumber;j++)
        {
            for(int k=j*sampleNumber;k<(j+1)*sampleNumber;k++)
            {
                signal[k][0]=signalTemp[k-j*sampleNumber][0];
                signal[k][1]=signalTemp[k-j*sampleNumber][1];
            }
        }

        for(int j=0;j<totalNumber;j++)
        {
            freqMove[j][0]=cos(2*PI*targetFd[i]*j/sampleRate);
            freqMove[j][1]=sin(2*PI*targetFd[i]*j/sampleRate);

            float tmp1=signal[j][0]*freqMove[j][0]-signal[j][1]*freqMove[j][1];
            float tmp2=signal[j][0]*freqMove[j][1]+signal[j][1]*freqMove[j][0];
            signal[j][0]=tmp1;
            signal[j][1]=tmp2;

            signalAll[j][0]=signalAll[j][0]+signal[j][0]; 
            signalAll[j][1]=signalAll[j][1]+signal[j][1]; 
        }

    }


    normrnd(0,pow(10.0,noisePower/10.0),systemNoise,totalNumber);

    // for(int i=0;i<totalNumber;i++)
    // {
    //     signalAll[i][0]=signalAll[i][0]+systemNoise[i][0];
    //     signalAll[i][1]=signalAll[i][1]+systemNoise[i][1];
    // }

    for(int i=0;i<pulseNumber;i++)
    {
        for(int j=0;j<blindNumber;j++)
            signalAll[i*sampleNumber+j][0]=signalAll[i*sampleNumber+j][1]=0.0;
    }

    _end=cpuSecond();
    double timeCost=_end-start;
	printf("Create Object Time Cost: %.6lfs\n",timeCost);

    fftw_free(signalTemp);
    fftw_free(signal);
    fftw_free(freqMove);
    fftw_free(delayNumber);
    fftw_free(targetFd);
    fftw_free(systemNoise);

    if(IO_mode==1)
    {
        FILE *fp=fopen("signal.dat", "w");
        for (int i=0;i<totalNumber;i++)
            fprintf(fp,"%f\n%f\n",signalAll[i][0],signalAll[i][1]);
        fclose(fp);
    }
}