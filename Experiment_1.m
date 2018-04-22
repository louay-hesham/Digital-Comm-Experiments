x = randi([0 1],1,(10^6));  
SNR = 0:2:30;
Average_power = ((1/10^6) * sum(x.^2)) % Average power
for i=1:length(SNR)         % Looping on the SNR values
    Rx_sequence=awgn( x, SNR(i),'measured' );    %Apply noise to the signal
    Rx_new=real(Rx_sequence)>0.5;       % Comparing signal+noise to threshhold to compute received signal
    [number,ratio] = biterr(x,Rx_new);  % Calculates number of errors and its ratio
    BER(i)= ratio;                     % Assign the probability of error of each SNR in matrix
end

semilogy(SNR,BER)   %Plot the BER curve against SNR
xlim([0 30])
title('BER VS SNR (simple detector)');
xlabel('SNR from 0dB to 30dB');
ylabel('BER');