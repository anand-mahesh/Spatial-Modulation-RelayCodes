clear all
            %this code can be generealised to any of the required mod
            %techniques given in the paper
Nt = 2;
Nr = 4;
Nu = 2;
M = 16;     %here working with QAM16
N = 4;

SNRdB = 0:1:15;
NUM = 10^4;
SER_SD = zeros(size(SNRdB));
SER = zeros(size(SNRdB));
count = zeros(size(SNRdB));
sigConstQAM16 = qammod((0:M-1).', M);
sigConstQ0 = qammod((0:N-1).', N);
sigConstQ1 = [1 -1 1i -1i];

spatialbits = 2;
eta = spatialbits + log2(M);

poss_ant_combQ = repmat([1 2], 2^(eta-1), 1);
poss_ant_combQAM16 = sort(repmat([1:2]', 2^(eta-2), 1));

poss_sig_symQAM16 = repmat((1:M).', 2^(spatialbits-1), 1);
poss_sig_symQ = repmat((1:N).', 2^(eta-2)/N, 2);
poss_sig_symQ(:, 1) = sort(poss_sig_symQ(:, 1));
poss_sig_symQ = repmat(poss_sig_symQ, 2, 1);

PossibleAntIndQAM16 = zeros(1, 2^(eta-1));
PossibleAntIndQ = zeros(Nu, 2^(eta-1));
PossibleAntIndQAM16(1, 1:2^(eta-1)) = sub2ind([Nu 2^(eta-1)], poss_ant_combQAM16.', 1:2^(eta-1));
for i = 1:Nu
    PossibleAntIndQ(i, 1:2^(eta-1)) =  sub2ind([Nu 2^eta], poss_ant_combQ(:, i).', 2^(eta-1)+1:2^eta);
end
% PossibleAntIndB = PossibleAntIndB + 2^(eta-1);

ESMconstdia = zeros(Nu, 2^eta);
ESMconstdia(PossibleAntIndQAM16(:, 1:2^(eta-1))) = sigConstQAM16(poss_sig_symQAM16);
for i = 1:Nu
    ESMconstdia(PossibleAntIndQ(i, 1:2^(eta-2))) = sigConstQ0(poss_sig_symQ(1:2^(eta-2), i));
    ESMconstdia(PossibleAntIndQ(i, 2^(eta-2)+1:2^(eta-1))) = sigConstQ1(poss_sig_symQ(2^(eta-2)+1:2^(eta-1), i));
end

for xx = 1:length(SNRdB)
    for i = 1:NUM
        pregendata = randi([0 (2^eta - 1)]); 
        x_t = ESMconstdia(:,pregendata+1);
        H = (1/sqrt(2))*(randn(Nr, Nt)+1j.*randn(Nr, Nt));

        snr = power(10,(SNRdB(xx)/10));
        sig_pwr = sum((abs(x_t)).^2)/max([mean(abs(sigConstQAM16)) mean(abs(sigConstQ0)) mean(abs(sigConstQ1))]);                                   % this is very important...
        noise_pwr = sig_pwr./snr;
        std_dev = power(noise_pwr,0.5);
        noise = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        
        y = H*x_t + noise;  %there is no spatial correlation
        
        %% sphere              
       [Q, R] = qr(H, 0);
       z = Q'*y;
       n = size(H,2);
       RADIUS = 20;
       TMPVAL        = zeros(n,1);
       SYMBSETSIZE   = length(ESMconstdia(1,:));
        global   RETVAL
        global SEARCHFLAG
        global SPHDEC_RADIUS
       for ii = 1:SYMBSETSIZE
          TMPVAL = ESMconstdia(:,ii);
          d=0;
          if ii==1
              SPHDEC_RADIUS= RADIUS;
          end
          if SPHDEC_RADIUS <= RADIUS
              RADIUS = SPHDEC_RADIUS;
          end
          sphdec_core(z, R,TMPVAL,n,d,RADIUS)

          if SEARCHFLAG > 0
             r = RETVAL;
          else
             r = zeros(2,1);
          end
       end
       if (r(1,1)~=x_t(1,1))|| (r(2,1)~=x_t(2,1))
           for check=1:length(ESMconstdia)
               if ESMconstdia(:,check)== r
                   break
               end    
           end
        siz2 = size(ESMconstdia,2);
        ML_Binary_Results2 = dec2bin(check-1,log2(siz2));
        BER_SMT_ML_SD =sum(dec2bin(pregendata,eta)~= ML_Binary_Results2)/eta;
        SER_SD(xx) = SER_SD(xx) +BER_SMT_ML_SD;
           count(xx) =count(xx)+1;
           
           siz = size(ESMconstdia, 2);
        [~,Idx_min_Error]= min(sum(abs((repmat(y,1,siz)-(H*ESMconstdia))).^2,1)) ;
        ML_Binary_Results = dec2bin(Idx_min_Error-1,log2(siz));
        BER_SMT_ML =(sum(dec2bin(pregendata,eta)~=ML_Binary_Results))/eta;
            
        SER(xx) = SER(xx) + BER_SMT_ML;
       end
    end
end
%%
SER_SD = SER_SD/NUM;
SER = SER/NUM;

figure
semilogy(SNRdB, SER_SD, 'x-', 'LineWidth', 1,'color','r', 'DisplayName','LOS');
hold on
semilogy(SNRdB,SER,'o-', 'LineWidth', 1,'color','g', 'DisplayName', 'Relay');
grid on 
xlabel('Es/No, dB') 
ylabel('Symbol Error Rate') 
xlabel('Average Eb/No,dB');
title('SER vs SNR of ESM-1 (Rayleigh Channel)');
axis([0 15 10^(-7) 1]);

 %%
 function sphdec_core(z, R,TMPVAL,layer, d,RADIUS)
 global   RETVAL
 global SEARCHFLAG
global SPHDEC_RADIUS
if (layer == 1)
        d = abs(z(1) - R(1,:)*TMPVAL)^2 + d;
        if (d <= RADIUS)
            RETVAL=TMPVAL;
            SPHDEC_RADIUS=d;
            SEARCHFLAG = 1;
            
        end
  else
        d = abs(z(layer) - R(layer,:)*TMPVAL)^2 + d;
        if d <=RADIUS
            sphdec_core(z, R,TMPVAL,layer-1,d,RADIUS);
        end
end

 end