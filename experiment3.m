
clc;
clear all;
close all;
% genetation of binary data vector
DataVector=randi([0 1],1,(10^3)); % generating 1000 because it takes too long lw 1^6 bta5d 45 d2ee2a
%

% Binary ASK demodulation (ASK modulation makes no change in the message zay ma hya  )
message=DataVector;
power_ask = mean(DataVector.^2);
SNR = 0:2:30 ;         %SNR with 2  steps
for i =1:1:length(SNR)
% Complex noise generation:
noise = sqrt( power_ask / (2 * db2pow(SNR(i))) ) * ( randn(size(DataVector)) + 1j *randn(size(DataVector)) );
reciv = message+noise ;    %adding the generated noise to message signal

 y=reciv>0.5;   % threshold b 0.5
% a comparision between the original bits with the detected ones
    [number,ratio] = biterr(message, y) ;  % calculating bit error rate
    BER_ASK(i)=[ratio];      %probability of error
end

%% PSK modulation

bit=[];
%if the transmitted bit is �1�, �1� is sent
% if it is a �0�,  �-1� is sent

for n=1:1:length(DataVector)
    if DataVector(n)==1;
       bit(n)=1;
    else
        bit(n)= -1;
    end
end

M=2;
message2 =  pskmod(double(DataVector) , M); % built in psk fl matlab

power_psk = mean(bit.^2);   % power of the transmitted signal

% Binary PSK demodulation
for i = 1:1:length(SNR)
%complex noise
noise = sqrt( power_psk / (2 *db2pow(SNR(i))) ) * ( randn(size(bit)) + 1j *randn(size(bit)) );
Rx_sequence = bit+noise ;     %adding the noise to the signal  

    modulated = message2+noise;  %adding the noise to the signal message

  bit1= Rx_sequence>=0;   %comparing with threshold =0

      dmodulated= pskdemod(modulated,2);           %demodulation by the built-in function � pskdemod �

[number,ratio] = biterr(DataVector,bit1);
    BER_PSK_manual(i)=[ratio]; % probability of error of each SNR in matrix
   [number,ratio1] = biterr(DataVector,dmodulated);
    BER1_PSK_function(i)=[ratio1]; % probability of error of each SNR in matrix
end

%%

% FSK

message2 = fskmod(DataVector,4,0.5,3,1.8); % built-in function  fskmod fl matlab.
snr=0:2:30;
%if the transmitted bit is �0�, �1� is sent
% if it is a �1�, the complex number �i� is sent

for i=1:length(DataVector)
   if DataVector(i)==1
      message(i)=1i;
   else
       message(i)=1;
   end
end
%getting power
power_fsk = mean(abs(message.^2));
power_fsk2 = mean(abs(message2.^2));

for j=1:length(snr)
                      %generating complex noise
noise=  sqrt(power_fsk/(db2pow(snr(j))))*(1/sqrt(2))*(randn(1,length(message)) +1i*randn(1,length(message)));
noise1=  (sqrt(power_fsk2/db2pow(snr(j))))*(1/sqrt(2))*(randn(1,length(message2)) +1i*randn(1,length(message2)));

                      %adding noise
    recieved = message+noise ;
    modulated=message2+noise1;
                      %demodulation by fskdemod function
    dmodulated2=fskdemod(modulated,4,0.5,3,1.8);
                      %FSK demodulation

decition= real(recieved)<= imag((recieved));     %if the real part of the signal equal or bigger than the imaginary part then the bit is zero and vice versa

[number,ratio] = biterr(DataVector,decition);
  BER_FSK_manual(j)=[ratio];% probability of error of each SNR in matrix
  [number,ratio1] = biterr(DataVector,dmodulated2);
    BER1_FSK_function(j)=[ratio1]; %probability of error of each SNR in matrix
end

%% plots
SNR=0:2:30;
 semilogy(SNR,[BER_ASK' BER_PSK_manual' BER_FSK_manual'])       %plot of manual codes
 grid on
 legend('ASK','PSK','FSK')
xlabel('SNR')
ylabel('BER')
title('BER vs SNR for manual codes')
figure
semilogy(SNR,[BER_ASK' BER1_PSK_function' BER1_FSK_function'])   %plot of function codes
grid on
 legend('ASK','PSK','FSK')
xlabel('SNR')
ylabel('BER')
title('BER vs SNR for function codes')


%%16-QAM
DataVector1= bi2de(reshape(DataVector,[],4));              %convert it in decimal
modulated= qammod(x1,16,0);        %16-QAM modulation

snr1=0:15:30;
for i=1:length(snr1)
     noisy_qam =  awgn( modulated, snr1(i),'measured' );  %applying noise to the signal

  noisy_plot=scatterplot(noisy_qam)   ;                     %%Plot the constellation diagram of the 16QAM at SNR values of 0 and 15 and 30 after adding  noise
     hold on
     scatterplot(modulated,1,0,'k*',noisy_plot)              %%Plot the original constellation diagram of the 16QAM at SNR =0 and 15 before adding noise
     title(['16-QAM, Binary Symbol Mapping',num2str(snr1(i))])


     demodulated=qamdemod(noisy_qam,16,0);                %demodulating the signal

      [number,ratio1] = biterr(DataVector1, demodulated) ;
    BER16(i)=[ratio1];%Assign the probability of error of each SNR in matrix
end

SNR=0:2:30;
for i=1:length(SNR)
     noisy_qam =  awgn( modulated, SNR(i),'measured' );  %applying noise to  the signal
      demodulated=qamdemod(noisy_qam,16,0);                %demodulating the signal

[No_errors,BER16(i)]=biterr(demodulated,DataVector1);
end
figure
semilogy(SNR,BER16)        %Plot the BER curve against SNR
xlabel(' snr1');
ylabel('BER1');
title('BER for 16-QAM')
