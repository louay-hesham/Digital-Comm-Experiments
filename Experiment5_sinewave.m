clc;
clear all;
close all;
%Sine wave Amplitude,Sampling frequency and sampling time.
A=2;
Fs=200e3;
Ts=1/Fs;
%x-coordinate (x-axis) from(0 till 5*10^-3)
t=0:Ts:5e-3-Ts;
%modulating signal
x= A*sin(2*pi*500*t);
len = length(t);
%plot the modulating signal
subplot(2,2,1);
plot(t,x,'r');
title('Modulating signal');
%Specify the length of the stair fn
delta = 0.2;
xn=0;
%start modulation
for i =1:len-1;
    if x(i)>xn(i)
        d(i)=1;
        xn(i+1)=xn(i)+delta;
    else
        d(i) =0;
        xn(i+1)=xn(i)-delta;
    end
end


%plot the signal after DM.
subplot(2,2,2);
stairs(t,xn);
title('signal after Delta modulation');
 %Demodulation
for i=1:d
    if d(i)>xn(i)
        d(i)=0;
        xn(i+1)=xn(i+1)-delta;
    else
        d(i)=1;
        xn(i+1)=xn(i)+delta;
    end
end
%plot demodulated signal
subplot(2,2,3);
plot(t,xn,'c');
title('signal after demodulation');

%BONUS
%Variable slope modulation
threshold = 15; %should be changed though
differential = diff(x);
vxn = 0;
for i =1:len-1;
    change = delta;
    if differential(i) > threshold
        change = 2 * delta;
    end
    if x(i)>vxn(i)
        vd(i)=1;
        vxn(i+1)=vxn(i)+change;
    else
        vd(i) =0;
        vxn(i+1)=vxn(i)-change;
    end
end
%Demodulation
for i=1:d
    if vd(i)>vxn(i)
        vd(i)=0;
        vxn(i+1)=vxn(i+1)-delta;
    else
        vd(i)=1;
        vxn(i+1)=vxn(i)+delta;
    end
end
%plot demodulated signal
subplot(2,2,4);
plot(t,vxn,'c');
title('signal after demodulation with variable delta');
%Demodulation using DPCM
xndpcm(1) = xn(1)
xndpcm(2) = xn(2)
for i = 3: len;
    xndpcm(i+1) = 2 * xndpcm(i) - xndpcm(i-1)


%Adjust the low pass filter
cut_off=1.5e3/Fs/2;
order=32;
h=fir1(order,cut_off);
%pass the demodulated signa; through the filter in the time domain.
con=conv(xn,h);
%plot the resulting signal after smoothing.
subplot(2,2,4);
plot(con);
title('Signal after smoothing');
%calc square error

