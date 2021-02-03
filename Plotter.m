close all

SNRdB = [0:1:20];

semilogy(SNRdB, ninetx, 'ro-', 'LineWidth', 1, 'DisplayName', '9 transmit');       
hold on
semilogy(SNRdB, fourtx, 'g*-','DisplayName', '4 transmit');
%semilogy(SNRdB, inter9, 'bo-', 'LineWidth', 1, 'DisplayName', 'Correlated Channel interleaved');
grid on  
ylabel('ABER') 
xlabel('Eb/N0,dB')
title('ABER vs SNR with Single Relay');
axis([0 15 10^(-5) 1]);