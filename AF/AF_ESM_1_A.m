clear all;
close all;

SNRdB = [0:1:20];
NUM = 10^4;
Nt = 8;
Nr_relay = 8;
Nr = 4;
Nu = 2;
M = 16;
N = 8;     
SER_SD = zeros(size(SNRdB));
SER_RD = zeros(size(SNRdB));
rho = 0.9;
K = 0.5;

%%
%ESM GENERATOR
% sigConstP = qammod((0:M-1).', M);
% sigConstS = [2 -2 2j -2j 2+2j 2-2j -2-2j -2+2j];    %assuming n = 8
% spatialbits = 3;
% eta = log2(M) + log2(N) + spatialbits;
% 
% all_ant_comb = nchoosek(1:Nt, Nu);
% all_ant_comb([2 5], :) = [];
% all_ant_comb = repmat(all_ant_comb, 2,1);
% poss_ant_comb = all_ant_comb(reshape(repmat(1:2^spatialbits, M*N, 1),2^eta,1) ,:);
% 
% poss_sig_symP = sort(repmat((1:M)', N, 1));
% poss_sig_symS = repmat((1:N)', 2^(spatialbits), 1);
% 
% PossibleAntInd = zeros(Nu, 2^eta);
% 
% for i = 1:Nu
%     PossibleAntInd(i, :) = sub2ind([Nt 2^eta], poss_ant_comb(:,i).', 1:2^eta);
% end
% 
% ESMconsdia = zeros(Nt, 2^eta);
% ESMconsdia(PossibleAntInd(1, 1:2^(eta-1))) = repmat(sigConstP(poss_sig_symP), 1,2^(eta-1)/size(poss_sig_symP,1));
% ESMconsdia(PossibleAntInd(1, 2^(eta-1)+1:2^eta)) = repmat(sigConstS(poss_sig_symS), 1, 2^(eta-1)/size(poss_sig_symS,1));
% ESMconsdia(PossibleAntInd(2, 1:2^(eta-1))) = repmat(sigConstS(poss_sig_symS), 1, 2^(eta-1)/size(poss_sig_symS,1));
% ESMconsdia(PossibleAntInd(2, 2^(eta-1)+1:2^eta)) = repmat(sigConstP(poss_sig_symP), 1, 2^(eta-1)/size(poss_sig_symP,1));


%% GENERALIZED ESM GENERATOR

sigConstP = qammod((0:M-1).', M);
sigConstS = [2 -2 2j -2j 2+2j 2-2j -2-2j -2+2j];
spatialbits = floor(log2(nchoosek(Nt,Nu)))+1;
eta = log2(M) + log2(N) + spatialbits;

unsort_ant_comb = nchoosek(1:Nt, Nu);
unsort_ant_comb(:,3) = unsort_ant_comb(:,2)-unsort_ant_comb(:,1);
all_ant_comb = sortrows(unsort_ant_comb, 3);
diff = size(unsort_ant_comb,1) - 2^(spatialbits-1);
all_ant_comb(1:diff,:) = [];
all_ant_comb(:,3) = [];
all_ant_comb = repmat(all_ant_comb, 2,1);
poss_ant_comb = all_ant_comb(reshape(repmat(1:2^spatialbits, M*N, 1),2^eta,1) ,:);

poss_sig_symP = sort(repmat((1:M)', N, 1));
poss_sig_symS = repmat((1:N)', 2^(spatialbits), 1);

PossibleAntInd = zeros(Nu, 2^eta);

for i = 1:Nu
    PossibleAntInd(i, :) = sub2ind([Nt 2^eta], poss_ant_comb(:,i).', 1:2^eta);
end

ESMconsdia = zeros(Nt, 2^eta);
ESMconsdia(PossibleAntInd(1, 1:2^(eta-1))) = repmat(sigConstP(poss_sig_symP), 1,2^(eta-1)/size(poss_sig_symP,1));
ESMconsdia(PossibleAntInd(1, 2^(eta-1)+1:2^eta)) = repmat(sigConstS(poss_sig_symS), 1, 2^(eta-1)/size(poss_sig_symS,1));
ESMconsdia(PossibleAntInd(2, 1:2^(eta-1))) = repmat(sigConstS(poss_sig_symS), 1, 2^(eta-1)/size(poss_sig_symS,1));
ESMconsdia(PossibleAntInd(2, 2^(eta-1)+1:2^eta)) = repmat(sigConstP(poss_sig_symP), 1, 2^(eta-1)/size(poss_sig_symP,1));


%% MONTE CARLO SIMULATION

Rtx = zeros(Nt, Nt);
for i = 1:Nt
    for j = 1:Nt
        Rtx(i,j) = rho^abs(i-j);
    end
end

for snrcount=1:length(SNRdB)  
    for tries = 0:NUM
        x = randi([0 2^eta-1]);
        x_s = ESMconsdia(:,x+1);
        
        %RAYLEIGH FADING
        H_sd_uncorr = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
        H_sr_uncorr = (1/sqrt(2))*(randn(Nr_relay,Nt) + 1i.*randn(Nr_relay,Nt));
        H_rd_uncorr = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
        
        %RICIAN FADING
        H_sd_uncorr = sqrt(K/(1+K))+ ones(size(H_sd_uncorr))+ sqrt(1/(1+K))*H_sd_uncorr;
        H_sr_uncorr = sqrt(K/(1+K))+ ones(size(H_sr_uncorr))+ sqrt(1/(1+K))*H_sr_uncorr;
        H_rd_uncorr = sqrt(K/(1+K))+ ones(size(H_rd_uncorr))+ sqrt(1/(1+K))*H_rd_uncorr;
 
        H_sd =  H_sd_uncorr*Rtx^(1/2);
        H_sr =  H_sr_uncorr*Rtx^(1/2);
        H_rd =  H_rd_uncorr*Rtx^(1/2);
        
        snr = 10^(SNRdB(snrcount)/10); %SNR
        sigpwr = sum((abs(x_s)).^2)/max([mean(abs(sigConstP)) mean(abs(sigConstS))]);
        noisepwr = sigpwr./snr;
        std_dev = sqrt(noisepwr);

        
        noise_sd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        noise_sr = std_dev.*(randn(Nr_relay,1)+1j*randn(Nr_relay,1)).*0.707;
        noise_rd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;

        %transmission from source to relay and destination
        y_sd = H_sd*x_s + noise_sd; %Destination received the signal 'y_sd' from Source 
        y_sr = H_sr*x_s + noise_sr; % Relay received the signal 'y_sr' from Source 
        
        siz = size(ESMconsdia,2);
        [~, Idx_min_Error] = min(sum(abs(repmat(y_sd,1,siz)-H_sd*ESMconsdia).^2,1));
        ML_Binary_Results = dec2bin(Idx_min_Error-1,log2(siz));
        BER_SMT_ML_SD =sum(dec2bin(x,eta)~= ML_Binary_Results)/eta;
        SER_SD(snrcount) = SER_SD(snrcount) + BER_SMT_ML_SD;
        
        %beta = sqrt(POW_S)/((abs(H_sr))^2 * POW_S + POW_N); % amplification factor 
        alpha = 0.5;
        
        
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

%% PLOTTER
figure
semilogy(SNRdB, SER_SD, 'o-', 'LineWidth', 1,'color','r', 'DisplayName','LOS');
hold on
semilogy(SNRdB, SER_RD,'o-', 'LineWidth', 1,'color','g', 'DisplayName', 'Relay');
grid on 
xlabel('Es/No, dB') 
ylabel('Symbol Error Rate') 
xlabel('Average Eb/No,dB');
title('SER vs SNR of ESM-1 (Rayleigh Channel)');
axis([0 15 10^(-7) 1]);
