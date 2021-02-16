clear all;
close all;

Nt  =8;
Nr_relay = 8;
rho = 0;
BLOCK_SIZE = 1023;
SNRdB = 0:1:15;
NUM = BLOCK_SIZE*10;
K=3.13;
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


%% GENERALISED ESM ENCODER

sigConstP = qammod((0:1:M-1).', M);
sigConstQ = sigConstP([2 7 12 13]);
sigConstP([1 2 3 7 9 11 12 13]) = [];
sigConstS = [-2+2j -2-2j 2j -2j 2+2j 2-2j -2 2].';

spatialbits1 = floor(log2(nchoosek(Nt,Nu)))+1;       %ETA is generated to match ESM1
eta = log2(M) + log2(N) + spatialbits1;
subspatialbits = 2;
spatialbits = eta - subspatialbits - 2*log2(N);

all_ant_comb = nchoosek(1:Nt, Nu);
all_ant_comb(:,3) = 1:size(nchoosek(1:Nt, Nu));
start = 1;
inc = 2^(spatialbits-1);
L1 = repmat(sortrows(repmat(all_ant_comb(start:start+inc-1, :), N.^2, 1),3),2,1);
L1(:,3) = [];
start = start + inc;
L2 = repmat(sortrows(repmat(all_ant_comb(start:start+inc-1, :), N.^2, 1),3),2,1);
L2(:,3) = [];
start = start + inc;
L3 = repmat(sortrows(repmat(all_ant_comb(start:start+inc-1, :), N.^2, 1),3),2,1);
L3(:,3) = [];
start = start + inc;
L4 = repmat(sortrows(repmat(all_ant_comb(1:2^spatialbits, :), 2^5, 1),3),2,1);
L4(:,3) = [];

poss_sig_symP = repmat(sort(repmat((1:N).', N, 1)),2^(eta-subspatialbits-log2(N))/N,1);
poss_sig_symS = repmat((1:N).', 2^(eta-subspatialbits)/N,1);
poss_sig_symQ = repmat(sort(repmat((1:4).', N, 1)),2^(eta-subspatialbits-log2(N/2))/N,1);

PossibleAntInd = zeros(Nu, 2^eta);

var = 2^(eta-subspatialbits);
for i = 1:Nu
    PossibleAntInd(i, 1:var) = sub2ind([Nt 2^eta], L1(:,i).', 1:var);
    PossibleAntInd(i, var+1:2*var) = sub2ind([Nt 2^eta], L2(:,i).', var+1:2*var);
    PossibleAntInd(i, 2*var+1:3*var) = sub2ind([Nt 2^eta], L3(:,i).', 2*var+1:3*var);
    PossibleAntInd(i, 3*var+1:4*var) = sub2ind([Nt 2^eta], L4(:,i).', 3*var+1:4*var);
end

ESMconsdia = zeros(Nt, 2^eta);

ESMconsdia(PossibleAntInd(1, 1:0.5*var)) = sigConstP(poss_sig_symP(1:0.5*var)); %L1
ESMconsdia(PossibleAntInd(2, 1:0.5*var)) = sigConstS(poss_sig_symS(1:0.5*var));
ESMconsdia(PossibleAntInd(1, 0.5*var+1:var)) = sigConstS(poss_sig_symS(0.5*var+1:var));
ESMconsdia(PossibleAntInd(2, 0.5*var+1:var)) = sigConstP(poss_sig_symP(0.5*var+1:var));

ESMconsdia(PossibleAntInd(1, var+1:1.5*var)) = sigConstP(poss_sig_symP(1:0.5*var));%L2
ESMconsdia(PossibleAntInd(2, var+1:1.5*var)) = sigConstS(poss_sig_symS(1:0.5*var));
ESMconsdia(PossibleAntInd(1, 1.5*var+1:2*var)) = sigConstS(poss_sig_symS(0.5*var+1:var));
ESMconsdia(PossibleAntInd(2, 1.5*var+1:2*var)) = sigConstP(poss_sig_symP(0.5*var+1:var));

ESMconsdia(PossibleAntInd(1, 2*var+1:2.5*var)) = sigConstP(poss_sig_symP(1:0.5*var));%L3
ESMconsdia(PossibleAntInd(2, 2*var+1:2.5*var)) = sigConstS(poss_sig_symS(1:0.5*var));
ESMconsdia(PossibleAntInd(1, 2.5*var+1:3*var)) = sigConstS(poss_sig_symS(0.5*var+1:var));
ESMconsdia(PossibleAntInd(2, 2.5*var+1:3*var)) = sigConstP(poss_sig_symP(0.5*var+1:var));

ESMconsdia(PossibleAntInd(1, 3*var+1:3.5*var)) = sigConstQ(poss_sig_symQ(1:0.5*var)); %L4
ESMconsdia(PossibleAntInd(2, 3*var+1:3.5*var)) = sigConstS(poss_sig_symS(1:0.5*var));
ESMconsdia(PossibleAntInd(1, 3.5*var+1:4*var)) = sigConstS(poss_sig_symS(0.5*var+1:var));
ESMconsdia(PossibleAntInd(2, 3.5*var+1:4*var)) = sigConstQ(poss_sig_symQ(0.5*var+1:var));
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
        H_sd = sqrt(K/(1+K))+ ones(size(H_sd))+ sqrt(1/(1+K))*H_sd;
        H_sd = H_sd*Rtx^(1/2);
        H_sr = (1/sqrt(2))*(randn(Nr_relay,Nt) + 1i.*randn(Nr_relay,Nt));
        H_sr = sqrt(K/(1+K))+ ones(size(H_sr))+ sqrt(1/(1+K))*H_sr;
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
                    H_rd = sqrt(K/(1+K))+ ones(size(H_rd))+ sqrt(1/(1+K))*H_rd;
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