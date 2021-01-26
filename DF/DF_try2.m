clear all;
close all;

THRESH = 0.09;
BLOCK_SIZE = 63;
SNRdB = 0:1:15;
NUM = BLOCK_SIZE*1000;
Nt  =4;
Nr = 4;
Nu = 2;
M = 16;
N = 8;     
SER_SD = zeros(size(SNRdB));
SER_RD = zeros(size(SNRdB));
SER_check=0;
counter = 0;

buffer1 = zeros(BLOCK_SIZE,1024);%calculating LOS path sums for ML detector
buffer2 = zeros(BLOCK_SIZE,1);%storing Tx symbol indices
buffer3 = zeros(BLOCK_SIZE*Nt,Nr);%storing H_sr matrices
buffer4 = zeros(Nr,BLOCK_SIZE);%storing y_sr matrices
%%
%ESM GENERATOR
sigConstP = qammod((0:M-1).', M);
sigConstS = [2 -2 2j -2j 2+2j 2-2j -2-2j -2+2j];    %assuming n = 8
spatialbits = 3;
eta = log2(M) + log2(N) + spatialbits;

all_ant_comb = nchoosek(1:Nt, Nu);
all_ant_comb([2 5], :) = [];
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
for snrcount=1:length(SNRdB)
    for tries = 1:NUM
        Idx = randi([0 2^eta-1]);
        buffer2(counter+1,1)= Idx;
        x_s = ESMconsdia(:,Idx+1);
        
        H_sd = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
        H_sr = (1/sqrt(2))*(randn(Nr,Nt) + 1i.*randn(Nr,Nt));
        buffer3((4*counter+1:4*counter+4),:)= H_sr;
 
        
        snr = 10^(SNRdB(snrcount)/10); %SNR
        sigpwr = sum((abs(x_s)).^2)/max([mean(abs(sigConstP)) mean(abs(sigConstS))]);
        noisepwr = sigpwr./snr;
        std_dev = sqrt(noisepwr);

        
        noise_sd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        noise_sr = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        
        
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
                    noise_rd = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
                    y_rd = H_rd*buffer4(:,k+1) + noise_rd;
                    gamma= H_rd*H_rd'+eye(Nr);
                    %decode both at the Rx
                    siz2 = size(ESMconsdia,2);
                    [~, Idx_min_Error2] = min(sum(abs(gamma^-0*(repmat(y_rd,1,siz2)-H_rd*buffer3((4*k+1:4*k+4),:)*ESMconsdia)).^2,1)+buffer1(k+1,:));
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
%%

SER_SD = SER_SD/NUM;
SER_RD = SER_RD/NUM;
figure
semilogy(SNRdB, SER_SD, 'o-', 'LineWidth', 1,'color','r', 'DisplayName','LOS');
hold on
semilogy(SNRdB, SER_RD, 'o-', 'LineWidth', 1,'color','g', 'DisplayName','relay');
grid on 
xlabel('Es/No, dB') 
ylabel('Symbol Error Rate') 
xlabel('Average Eb/No,dB');
title('SER vs SNR of ESM-1 (Rayleigh Channel)');
axis([0 15 10^(-7) 1]);
