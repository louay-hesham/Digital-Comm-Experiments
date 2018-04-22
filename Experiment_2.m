x = randi([0 1],1,(10^6));  
SNR = 0:2:30;               
m = 10;
sampling_number = 10;
threshold=5*sqrt(1/10);
%x_waveform = nan(1,10^7);
%for i=1:length(x)
 %   for j=1:length(x)
  %      if (i ~= j)
   %         x_waveform(i) = [x(i) x(j)];
    %    end
   % end
% end
x_waveform=reshape(repmat(x,m,1),1,[]).*sqrt(1/10);  %reapeat matrix m 10 times then reshaping it to 1*10^7 vector
Average_power = ((1/10^7) * sum(x_waveform.^2)) % Average power of transimtter
for i=1:length(SNR)         % Looping on the SNR values
    Rx_sequence=awgn(x_waveform, SNR(i),'measured' );    %Apply noise to the signal
    Rx_sequence=reshape(Rx_sequence,10,[]);
    ten_1s = ones(sampling_number,1); % needed for the next line
    hMF = conv2(Rx_sequence,ten_1s,'valid'); % equivelant to hMF = c*(s1(T-t)-s2(T-t))
    ret_conv = hMF>=threshold; % Get 1's and 0's
    [number1,ratio1]=biterr(ret_conv,x);  % Calculates number of errors and its ratio
    BER_conv(i) = ratio1;
    g=ones(10,1); %g=s1
    gF=(Rx_sequence'*g)'; %multiplying received signal element by element and add resultant samples together
    ret_corr=gF>=threshold; % Get 1's and 0's
    [number2,ratio2] = biterr(ret_corr,x);  % Calculates number of errors and its ratio
    BER_cor(i)= ratio2; % Assign the probability of error of each SNR in matrix for correlator
end
% first Experiment
for i=1:length(SNR)         % Looping on the SNR values
    Rx_sequence=awgn( x, SNR(i),'measured' );    %Apply noise to the signal
    Rx_new=real(Rx_sequence)>0.5;       % Comparing signal+noise to threshhold to compute received signal
    [number,ratio] = biterr(x,Rx_new);  % Calculates number of errors and its ratio
    BER(i)= ratio;                     % Assign the probability of error of each SNR in matrix
end

%%% plots
semilogy(SNR,BER_conv)  %Plot the BER curve against SNR
xlim([0 30])
title('BER Of Matched Filter VS SNR');
xlabel('SNR from 0dB to 30dB');
ylabel('BER');

figure
semilogy(SNR,BER_cor)  %Plot the BER curve against SNR
xlim([0 30])
title('BER Of Correlator VS SNR');
xlabel('SNR from 0dB to 30dB');
ylabel('BER');
figure

semilogy(SNR,[BER(1:16)'  BER_conv'])%Plot the BER curve against SNR
legend('Simple detector','Matched filter and Correlator')
xlim([0 30])
title('Comparison between Matched filter and simple detector(BER VS SNR) ');
xlabel('SNR from 0dB to 30dB');
ylabel('BER')

