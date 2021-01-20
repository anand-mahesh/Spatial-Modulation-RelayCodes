clear all;
close all;

Nt = 4;
Nu = 2;
Nr = 4;

%rho = 0.9;
M = 16;
N = M/2;

SNRdB = [0:1:15];
NUM = 10e4;
SER_SD = zeros(size(SNRdB));
SER_RD = zeros(size(SNRdB));

sigConstP = qammod((0:1:M-1).', M);
sigConstQ = sigConstP([2 7 12 13]);
sigConstP([1 2 3 7 9 11 12 13]) = [];
sigConstS = [-2+2j -2-2j 2j -2j 2+2j 2-2j -2 2].';

subspatialbits = 2;
spatialbits = 2;

eta = subspatialbits + spatialbits + 2*log2(N);

all_ant_comb = nchoosek(1:Nt, Nu);
L1 = repmat(sort(repmat(all_ant_comb([1 6], :), N.^2, 1)),2,1);
L2 = repmat(sort(repmat(all_ant_comb([2 5], :), N.^2, 1)),2,1);
L3 = repmat([repmat(all_ant_comb(3, :) ,64,1); repmat(all_ant_comb(4, :) ,64,1)], 2, 1);          %this is the fucking error
combindex = repmat(reshape(repmat([1 6 2 5], 2^5,1), 4*2^5, 1),2,1);
L4 = all_ant_comb(combindex, :);

poss_sig_symP = repmat(sort(repmat((1:N).', N, 1)),4,1);
poss_sig_symS = repmat((1:N).', 32, 1);
poss_sig_symQ = repmat(sort(repmat((1:4).', N, 1)),8,1);

PossibleAntInd = zeros(Nu, 2^eta);

for i = 1:Nu
    PossibleAntInd(i, 1:256) = sub2ind([Nt 2^eta], L1(:,i).', 1:256);
    PossibleAntInd(i, 257:512) = sub2ind([Nt 2^eta], L2(:,i).', 257:512);
    PossibleAntInd(i, 513:768) = sub2ind([Nt 2^eta], L3(:,i).', 513:768);
    PossibleAntInd(i, 769:1024) = sub2ind([Nt 2^eta], L4(:,i).', 769:1024);
end

ESMconsdia = zeros(Nt, 2^eta);

ESMconsdia(PossibleAntInd(1, 1:128)) = sigConstP(poss_sig_symP(1:128)); %L1
ESMconsdia(PossibleAntInd(2, 1:128)) = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 129:256)) = sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 129:256)) = sigConstP(poss_sig_symP(129:256));

ESMconsdia(PossibleAntInd(1, 257:384)) = sigConstP(poss_sig_symP(1:128));%L2
ESMconsdia(PossibleAntInd(2, 257:384)) = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 385:512)) = sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 385:512)) = sigConstP(poss_sig_symP(129:256));

ESMconsdia(PossibleAntInd(1, 513:640)) = sigConstP(poss_sig_symP(1:128));%L3
ESMconsdia(PossibleAntInd(2, 513:640)) = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 641:768)) = sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 641:768)) = sigConstP(poss_sig_symP(129:256));

ESMconsdia(PossibleAntInd(1, 769:896)) = sigConstQ(poss_sig_symQ(1:128)); %L4
ESMconsdia(PossibleAntInd(2, 769:896)) = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 897:1024)) = sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 897:1024)) = sigConstQ(poss_sig_symQ(129:256));

%%
for snrcount=1:length(SNRdB)  
    for tries = 0:NUM
        x = randi([0 2^eta-1]);
        x_s = ESMconsdia(:,x+1);
        
        H_sd = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
        H_sr = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
        H_rd = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
 
        
        snr = 10^(SNRdB(snrcount)/10); %SNR
        sigpwr = sum((abs(x_s)).^2)/max([mean(abs(sigConstP)) mean(abs(sigConstS)) mean(abs(sigConstQ))]);
        noisepwr = sigpwr./snr;
        std_dev = sqrt(noisepwr);

        
        noise_sd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        noise_sr = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        noise_rd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;

        %transmission from source to relay and destination
        y_sd = H_sd*x_s + noise_sd; %Destination received the signal 'y_sd' from Source 
        y_sr = H_sr*x_s + noise_sr; % Relay received the signal 'y_sr' from Source 
        
        siz = size(ESMconsdia,2);
        [~, Idx_min_Error] = min(sum(abs(repmat(y_sd,1,siz)-H_sd*ESMconsdia).^2,1));
        ML_Binary_Results = dec2bin(Idx_min_Error-1,log2(siz));
        BER_SMT_ML_SD =sum(dec2bin(x,eta)~= ML_Binary_Results)/eta;
        SER_SD(snrcount) = SER_SD(snrcount) + BER_SMT_ML_SD;
        
        alpha = 0.5;        %AMPLIFICATION FACTOR
        
        y_rd = alpha.*H_rd*y_sr + noise_rd;
        gamma= alpha^2.*H_rd*H_rd'+eye(Nr);
        [~, Idx_min_Error1] = min(sum(abs(repmat(y_sd,1,siz)-H_sd*ESMconsdia).^2,1)+sum((abs(gamma^-0.5*(repmat(y_rd,1,siz)-(alpha*H_rd*H_sr*ESMconsdia)))).^2,1));
        ML_Binary_Results1 = dec2bin(Idx_min_Error1-1,log2(siz));
        BER_SMT_ML_RD = sum(dec2bin(x,eta)~= ML_Binary_Results1)/eta;
        SER_RD(snrcount) = SER_RD(snrcount) + BER_SMT_ML_RD;        
    end
end
SER_SD = SER_SD/NUM;
SER_RD = SER_RD/NUM;

%%
figure
semilogy(SNRdB, SER_SD, 'o-', 'LineWidth', 1,'color','r', 'DisplayName','LOS');
hold on
semilogy(SNRdB, SER_RD,'o-', 'LineWidth', 1,'color','g', 'DisplayName', 'Relay');
grid on 
xlabel('Es/No, dB') 
ylabel('Symbol Error Rate') 
xlabel('Average Eb/No,dB');
title('SER vs SNR of ESM-2 (Rayleigh Channel)');
axis([0 15 10^(-7) 1]);