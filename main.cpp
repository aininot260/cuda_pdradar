#include <cstdio>
#include <fftw3.h>

#include "utils.h"
#include "pd_radar.h"

int main()
{
    float sampleRate=2.0e6;  //采样频率
    float timeWidth=4.0e-5; //信号时宽
    float bandWidth=1.0e6;  //信号带宽
    float f0=0.0;  //中心频率
    float PRT=4.096e-3;  //雷达发射脉冲重复周期,240us对应1/2*240*300=36000米
    float wavegate=4.096e-3;  //波门大小
    int pulseNumber=32;  //回波脉冲数
    int targetNumber=3;  //目标个数
    float RF=3.140e9/2;  //雷达射频 电磁波频率
    int noisePower=-12;  //噪声功率,目标功率0dB
    bool pcFlag=0;  //0表示频域脉压,1表示时域脉压,时域脉压不加窗
    double parameters1=0.0;  //PC加窗参数
    int windowType1=4;  //PC加窗类型
    double parameters2=0.0;  //MTD加窗参数
    int windowType2=4;  //MTD加窗类型
    int protectNumber=1;  //保护单元的个数
	int refWindowSize=16;  //CFAR窗长
    double alpha=0.0;  //LogSquare参数
    double pf=5.0e-4;  //CFAR阈值
    int cfarType=1;  //CFAR类型
	int detectorType=4;  //CFAR探测类型
    int blindNumber=set_blindNumber(sampleRate,timeWidth);  //设置遮挡点的个数
    int sampleNumber=set_sampleNumber(sampleRate,PRT);  //设置一个脉冲的采样点
    int totalNumber=set_totalNumber(sampleNumber,pulseNumber);  //设置总的采样点

    fftwf_complex *coeff,*signalAll;
    float *cfar,*threshold;
    malloc_IO_memory(&signalAll,&coeff,&cfar,&threshold);
    int *signalPower=(int *)fftw_malloc(targetNumber*sizeof(int));  //目标功率,无量纲
    int *targetDistance=(int *)fftw_malloc(targetNumber*sizeof(int));  //目标距离,单位m
    int *targetVelocity=(int *)fftw_malloc(targetNumber*sizeof(int));  //目标径向速度 单位m/s
    fftwf_complex *chirp=(fftwf_complex *)fftw_malloc(blindNumber*sizeof(fftwf_complex));
    fftwf_complex *pc=(fftwf_complex *)fftw_malloc(totalNumber*sizeof(fftwf_complex));
    fftwf_complex *mti=(fftwf_complex *)fftw_malloc((totalNumber-sampleNumber)*sizeof(fftwf_complex));
    fftwf_complex *mtd=(fftwf_complex *)fftw_malloc((totalNumber)*sizeof(fftwf_complex));
    float *freqWeight=(float *)fftw_malloc((totalNumber)*sizeof(float));

    targetVelocity[0]=80,targetVelocity[1]=160,targetVelocity[2]=240;
    signalPower[0]=signalPower[1]=signalPower[2]=1;
    targetDistance[0]=50000,targetDistance[1]=80000,targetDistance[2]=250000;
    for(int i=0;i<totalNumber;i++) freqWeight[i]=1.0;

    set_gpu_acceletate(1);
    set_IO_mode(0);  //0:循环模式 1:文件输出 2：全部同步到内存
    set_stream_number(1);  //设置GPU程序使用的流，该设置不能超过最大流总数
    initialize(parameters1,parameters2,windowType1,windowType2,freqWeight,sampleRate,wavegate);
    
    call_LFM(chirp,coeff,timeWidth,bandWidth,f0);
    generate_object(chirp,signalAll,targetNumber,signalPower,targetDistance,targetVelocity,RF,noisePower);
    
    clean_timer();

    for(int i=0;i<1000;i++)
    {
        call_PC(signalAll,coeff,pc,pcFlag);
        call_MTI(pc,mti);
        call_MTD(mti,mtd,freqWeight);
        call_CFAR(mtd,cfar,threshold,protectNumber,refWindowSize,alpha,pf,cfarType,detectorType);
    }

    sync_time();
    get_total_time();

    fftw_free(pc);
    fftw_free(mti);
    fftw_free(mtd);
    fftw_free(freqWeight);
    fftw_free(chirp);
    fftw_free(signalPower);
    fftw_free(targetDistance);
    fftw_free(targetVelocity);

    free_IO_memory(signalAll,coeff,cfar,threshold);
    release();

    return 0;
}
