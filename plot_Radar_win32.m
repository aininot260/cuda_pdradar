clear all;
SampleRate=2.0e6;
BandWidth=1.0e6;
PulseNumber=32;
TimeWidth=4.0e-5;
PRT=4.096e-3;
SampleNumber=PRT*SampleRate;
BlindNumber=SampleRate*TimeWidth;
t=-TimeWidth/2:1/SampleRate:TimeWidth/2-1/SampleRate;

load('./build/Release/LFM.dat');
LFM_c=LFM(1:2:end) + LFM(2:2:end)*1i;
figure(10);
subplot(311);
plot(t,real(LFM_c));title('LFM信号实部波形');xlabel('t/s');ylabel('幅度');
subplot(312);
plot(t,imag(LFM_c));title('LFM信号虚部波形');xlabel('t/s');ylabel('幅度');
f=linspace(-SampleRate/2,SampleRate/2,BlindNumber);
LFM_fft=fftshift(fft(LFM_c));
subplot(313);
plot(f,abs(LFM_fft));title('LFM信号频谱图');xlabel('f/Hz');ylabel('幅度');

load('./build/Release/signal.dat');
signal_c=signal(1:2:end) + signal(2:2:end)*1i;
figure(9);
plot(abs(signal_c));title('回波信号(展开)');

load('./build/Release/pc.dat');
pc_c=pc(1:2:end) + pc(2:2:end)*1i;
figure(7);
plot(abs(pc_c));title('PC 结果(展开)');

for i=1:PulseNumber
      pc_c_3d(i,1:SampleNumber)=pc_c((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(8);
% subplot(121);
mesh(abs(pc_c_3d));title('PC 结果(立体)');xlabel('距离单元');ylabel('多普勒通道');zlabel('幅度');

PulseNumber=PulseNumber-1;

load('./build/Release/mti.dat');
mti_c=mti(1:2:end) + mti(2:2:end)*1i;
%figure(5);plot(abs(mti_c));title('MTI 结果(展开)');

for i=1:PulseNumber
      mti_c_3d(i,1:SampleNumber)=mti_c((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(6);
subplot(121);
mesh(abs(mti_c_3d));title('MTI 结果(立体)');xlabel('距离单元');ylabel('多普勒通道');zlabel('幅度');

load('./build/Release/mtd.dat');
mtd_c=mtd(1:2:end) + mtd(2:2:end)*1i;
%figure(1);plot(abs(mtd_c));title('MTD 结果(展开)');

for i=1:PulseNumber
      mtd_c_3d(i,1:SampleNumber)=mtd_c((i-1)*SampleNumber+1:i*SampleNumber);
end
% figure(2);
subplot(122);
mesh(abs(mtd_c_3d));title('MTD 结果(立体)');xlabel('距离单元');ylabel('多普勒通道');zlabel('幅度');

load('./build/Release/cfar.dat');
%figure(3);plot(abs(cfar));title('CFAR 结果(展开)');
for i=1:PulseNumber
      cfar_3d(i,1:SampleNumber)=cfar((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(4);
mesh(abs(cfar_3d));title('CFAR 结果(立体)');xlabel('距离单元');ylabel('多普勒通道');zlabel('幅度');

load('./build/Release/threshold.dat');
figure(11);plot(abs(cfar));hold on;plot(abs(threshold));title('CFAR 门限');
