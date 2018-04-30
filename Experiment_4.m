%Experiment 4
%PCM modulation
clc
clear 
%Generating sinusoidal wave
a=1; %Amplitude
f=2; %frequency
fs=4000; %sampling frequency
ts=0:1/fs:1;
no_o_bits = 'number of bits =  ';
n = input(no_o_bits);
m=2*n+1;
%number of levels
levels=2*n; 
delta= (2*a)/levels;
xs=a*sin(2*pi*f*ts);
%Quantization
%a zero vector of size 4001, the same lenght of ts
quan= zeros(1,4001);
levels_bits = zeros(1,1);
bits = zeros(1,1);
%quantization using fi function which return a fixed point object
%colon operator reshapes all elements into a single column vector.
quan(:) =double(fi(xs,1,m,n));
%ploting
figure(1)
subplot(2,1,1);
plot(ts,xs,'g')
title('Input Signal')

subplot(2,1,2); 
plot(ts,quan(:),'r')
title('Quantized Signal')

%Encoding
% saving unique number of quantized values in variable z for mapping
z= unique(quan);
for i = 1:length(z)
        % Encoding quantization level 
        levels_bits(i) = abs(z(i)*(2^n)+(2^n));
        %disp(levels_bits(i));
        % Converting the output integer values to binary using de2bi
        bits = de2bi(levels_bits(i),n+2,'left-msb');
        %disp(bits);      
end

%Calculating MSE
MSE = 0;
for i = 1:4001
        MSE = MSE+((quan(i)-xs(i))^2);
end
 MSE = (1/4001)*MSE;
 % Printing out the MSE
 fprintf('The value of MSE when n = %d is %f\n',n,MSE);

%==================================== Part two ====================================

% reconstruction from oversampling
figure
t=0:0.001:1;                    % time signal
y=2*cos(2*pi*5*t);
[B,A] = butter(3,1000/100000,'low');        % butter fly filter
zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
    zero_added_signal(i*10)=y(i);
end
zero_added_signal(1:9)=[];
% Adding zeros enhances the signal display and don't change the
% spectrum, it changes sampling freq. only
t=linspace(0,1,length(zero_added_signal));
%filter --> filters a data sequance using a digital filter 
%filters data with finite impulse response
filtered_signal = filter(B,A,zero_added_signal);
figure(2);
plot (t,filtered_signal,'r');
xlabel('time')
ylabel('oversampled signal')

% construction from minimum sampling
%fm = 5 for minimun sampling fs = 2fm --> fs =10 ts =1/fs 
figure
t=0:0.1:1;   % -------------------------------------------------------------------------> here 
y=2*cos(2*pi*5*t);
[B,A] = butter(10,0.1,'low');        
zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
    zero_added_signal(i*10)=y(i);
end
zero_added_signal(1:9)=[];
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
figure(3)
plot(t,filtered_signal,'r')
xlabel('time')
ylabel('minimum sampled signals')
s=fft(filtered_signal);
s=fftshift(s);
fs=100;     % why 100? ----------------------------------------------------------------> something wrong according to Eng esraa >
freq=linspace(-fs/2,fs/2,length(s));
figure(4)
plot (freq,abs(s));
xlabel('freq')
ylabel('magnitude of minimum sampled signal')

% construction from undesampling sampling
figure
t=0:0.2:1;
y=2*cos(2*pi*5*t);
[B,A] = butter(10,0.2,'low');
zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
    zero_added_signal(i*10)=y(i);
end
zero_added_signal(1:9)=[];
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
figure(5);
plot(t,filtered_signal,'r')
xlabel('time')
ylabel('undersampled signals')