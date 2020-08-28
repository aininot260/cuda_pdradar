clear all;
SampleRate=2.0e6;
BandWidth=1.0e6;
PulseNumber=32;
TimeWidth=4.0e-5;
PRT=4.096e-3;
SampleNumber=PRT*SampleRate;
BlindNumber=SampleRate*TimeWidth;
t=-TimeWidth/2:1/SampleRate:TimeWidth/2-1/SampleRate;


load('./build/LFM.dat');
LFM_c=LFM(1:2:end) + LFM(2:2:end)*1i;
figure(10);
subplot(311);
plot(t,real(LFM_c));title('LFM�ź�ʵ������');xlabel('t/s');ylabel('����');
subplot(312);
plot(t,imag(LFM_c));title('LFM�ź��鲿����');xlabel('t/s');ylabel('����');
f=linspace(-SampleRate/2,SampleRate/2,BlindNumber);
LFM_fft=fftshift(fft(LFM_c));
subplot(313);
plot(f,abs(LFM_fft));title('LFM�ź�Ƶ��ͼ');xlabel('f/Hz');ylabel('����');

load('./build/signal.dat');
signal_c=signal(1:2:end) + signal(2:2:end)*1i;
figure(9);
plot(abs(signal_c));title('�ز��ź�(չ��)');

load('./build/pc.dat');
pc_c=pc(1:2:end) + pc(2:2:end)*1i;
figure(7);
plot(abs(pc_c));title('PC ���(չ��)');

for i=1:PulseNumber
      pc_c_3d(i,1:SampleNumber)=pc_c((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(8);
mesh(abs(pc_c_3d));title('PC ���(����)');

PulseNumber=PulseNumber-1;

load('./build/mti.dat');
mti_c=mti(1:2:end) + mti(2:2:end)*1i;
%figure(5);plot(abs(mti_c));title('MTI ���(չ��)');

for i=1:PulseNumber
      mti_c_3d(i,1:SampleNumber)=mti_c((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(6);
mesh(abs(mti_c_3d));title('MTI ���(����)');

load('./build/mtd.dat');
mtd_c=mtd(1:2:end) + mtd(2:2:end)*1i;
%figure(1);plot(abs(mtd_c));title('MTD ���(չ��)');

for i=1:PulseNumber
      mtd_c_3d(i,1:SampleNumber)=mtd_c((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(2);
mesh(abs(mtd_c_3d));title('MTD ���(����)');

load('./build/cfar.dat');
%figure(3);plot(abs(cfar));title('CFAR ���(չ��)');
for i=1:PulseNumber
      cfar_3d(i,1:SampleNumber)=cfar((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(4);
mesh(abs(cfar_3d));title('CFAR ���(����)');

load('./build/threshold.dat');
figure(11);plot(abs(cfar));hold on;plot(abs(threshold));title('CFAR ����');
