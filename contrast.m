clear all;
SampleRate=2.0e6;
BandWidth=1.0e6;
PulseNumber=32;
TimeWidth=4.0e-5;
PRT=4.096e-3;
SampleNumber=PRT*SampleRate;
BlindNumber=SampleRate*TimeWidth;
t=-TimeWidth/2:1/SampleRate:TimeWidth/2-1/SampleRate;

load('./cpu_results/pc.dat');
pc_c=pc(1:2:end) + pc(2:2:end)*1i;
load('./gpu_results/pc.dat');
pc_d=pc(1:2:end) + pc(2:2:end)*1i;
figure(1);
plot(abs(pc_c-pc_d));title('PC 结果对比');

load('./cpu_results/mti.dat');
mti_c=mti(1:2:end) + mti(2:2:end)*1i;
figure(21);
plot(abs(mti_c));title('MTIs CPU结果');
load('./gpu_results/mti.dat');
mti_d=mti(1:2:end) + mti(2:2:end)*1i;
figure(2);
plot(abs(mti_c-mti_d));title('MTI 结果对比');

load('./cpu_results/mtd.dat');
mtd_c=mtd(1:2:end) + mtd(2:2:end)*1i;
figure(31);
plot(abs(mtd_c));title('MTD CPU结果');
load('./gpu_results/mtd.dat');
mtd_d=mtd(1:2:end) + mtd(2:2:end)*1i;
figure(32);
plot(abs(mtd_d));title('MTD GPU结果');
figure(3);
plot(abs(mtd_c-mtd_d));title('MTD 结果对比');

load('./cpu_results/cfar.dat');
cfar_c=cfar;
load('./gpu_results/cfar.dat');
cfar_d=cfar;
figure(4);
plot(abs(cfar_c-cfar_d));title('CFAR 结果对比');