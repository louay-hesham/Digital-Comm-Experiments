clc;
clear all;
close all;
%Sine wave Amplitude,Sampling frequency and sampling time.
A=2;
Fs=10000;
Ts=1/Fs;
%x-coordinate (x-axis) from(0 till 5*10^-3)
t=0:Ts:0.1;
%modulating signal
x= A*sin(2*pi*Fs/500*t);
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
 %Demodualtion 
for i=1:d
    if d(i)>xn(i)
        d(i)=0;
        xn(i+1)=xn(i)-delta;
    else
        d(i)=1;
        xn(i+1)=xn(i)+delta;
    end
end
%plot demodulated signal
subplot(2,2,3);
plot(t,xn,'c');
title('signal after demodulation');
%Adjust the low pass filter
order=3;
frame=81;
h= sgolayfilt(xn,order,frame);
%plot the resulting signal after smoothing.
subplot(2,2,4);
plot(t,h);
title('Signal after smoothing');
%calc square error
err = immse(x,h );
fprintf('\n The mean-squared error is %0.4f\n', err);