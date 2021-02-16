clear all;
close all;
Nt  =8;
Nr_relay = 8;
rho = 0;
BLOCK_SIZE = 1023;
SNRdB = 0:1:15;
NUM = BLOCK_SIZE*10;
%%
THRESH = 0.05;
Nr = 4;
Nu = 2;
M = 16;
N = 8;     
SER_SD = zeros(size(SNRdB));
SER_RD = zeros(size(SNRdB));
SER_check=0;
counter = 0;
rho = 0;

%%
%ESM generator

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
%%
buffer1 = zeros(BLOCK_SIZE,2^eta);%calculating LOS path sums for ML detector
buffer2 = zeros(BLOCK_SIZE,1);%storing Tx symbol indices
buffer3 = zeros(BLOCK_SIZE*Nt,Nr_relay);%storing H_sr matrices
buffer4 = zeros(Nr_relay,BLOCK_SIZE);%storing y_sr matrices

Rtx = zeros(Nt, Nt);
for i = 1:Nt
    for j = 1:Nt
        Rtx(i,j) = rho^abs(i-j);
    end
end
%%
for snrcount=1:length(SNRdB)
    for tries = 1:NUM
        Idx = randi([0 2^eta-1]);
        buffer2(counter+1,1)= Idx;
        x_s = ESMconsdia(:,Idx+1);
        
        H_sd = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
        H_sd = H_sd*Rtx^(1/2);
        H_sr = (1/sqrt(2))*(randn(Nr_relay,Nt) + 1i.*randn(Nr_relay,Nt));
        H_sr = H_sr*Rtx^(1/2);
        buffer3((Nt*counter+1:Nt*counter+Nt),:)= H_sr;
 
        
        snr = 10^(SNRdB(snrcount)/10); %SNR
        sigpwr = sum((abs(x_s)).^2)/max([mean(abs(sigConstP)) mean(abs(sigConstS))]);
        noisepwr = sigpwr./snr;
        std_dev = sqrt(noisepwr);

        
        noise_sd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        noise_sr = std_dev.*(randn(Nr_relay,1)+1j*randn(Nr_relay,1)).*0.707;
             
        if counter <= BLOCK_SIZE
            %transmission from source to destination
            y_sd = H_sd*x_s + noise_sd; %Destination received the signal 'y_sd' from Source 
        
            siz = size(ESMconsdia,2);
            buffer1(counter+1,:)=sum(abs(repmat(y_sd,1,siz)-H_sd*ESMconsdia).^2,1);
            [~, Idx_min_Error] = min(buffer1(counter+1,:));
            ML_Binary_Results = dec2bin(Idx_min_Error-1,log2(siz));
            BER_SMT_ML_SD =sum(dec2bin(Idx,eta)~= ML_Binary_Results)/eta;
            SER_SD(snrcount) = SER_SD(snrcount) +BER_SMT_ML_SD;
            
        
            %transmission from source to relay
            y_sr = H_sr*x_s + noise_sr; 
            buffer4(:,counter+1)=y_sr;
            siz1 = size(ESMconsdia,2);
            [~, Idx_min_Error1] = min(sum(abs(repmat(y_sr,1,siz1)-H_sr*ESMconsdia).^2,1));
            ML_Binary_Results1 = dec2bin(Idx_min_Error1-1,log2(siz1));
            BER_SMT_ML_SR =sum(dec2bin(Idx,eta)~= ML_Binary_Results1)/eta;
            SER_check =SER_check+BER_SMT_ML_SR;
        end
        if counter == BLOCK_SIZE
            if SER_check/BLOCK_SIZE <= THRESH
                for k =1:BLOCK_SIZE    
                    H_rd = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
                    H_rd = H_rd*Rtx^(1/2);
                    noise_rd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
                    y_rd = H_rd*buffer4(:,k+1) + noise_rd;
                    %decode both at the Rx
                    siz2 = size(ESMconsdia,2);
                    [~, Idx_min_Error2] = min(sum(abs((repmat(y_rd,1,siz2)-H_rd*buffer3((Nt*k+1:Nt*k+Nt),:)*ESMconsdia)).^2,1)+buffer1(k+1,:));
                    ML_Binary_Results2 = dec2bin(Idx_min_Error2-1,log2(siz2));
                    BER_SMT_ML_RD =sum(dec2bin(buffer2(k+1,1),eta)~= ML_Binary_Results2)/eta;
                    SER_RD(snrcount) = SER_RD(snrcount) +BER_SMT_ML_RD;
                end
            else 
               for k =1:BLOCK_SIZE
                   %decode only SD at the Rx
                    siz2 = size(ESMconsdia,2);
                    [~, Idx_min_Error2] = min(buffer1(k+1,:));
                    ML_Binary_Results2 = dec2bin(Idx_min_Error2-1,log2(siz2));
                    BER_SMT_ML_RD =sum(dec2bin(buffer2(k+1,1),eta)~= ML_Binary_Results2)/eta;
                    SER_RD(snrcount) = SER_RD(snrcount) +BER_SMT_ML_RD;
               
               end
            end
            counter=0;
            SER_check =0;
        end
     counter=counter+1;   
    end
end
%% PLOTTER
SER_SD = SER_SD/NUM;
SER_RD = SER_RD/NUM;
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
