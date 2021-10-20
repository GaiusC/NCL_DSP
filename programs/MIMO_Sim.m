%Author: J.Chen
%Matlab Version: R2021a
clear, clc, close all
tic
%===Reset random number===
RN=sum(100*clock);
RS=RandStream('mt19937ar','seed',RN);
RandStream.setGlobalStream(RS);

%===Simulation parameter definitions===
%---16QAM parameter---
M = 16;% M for 16 QAM
symbol_bits = log2(M);% 4 bits for a symbol in 16 QAM
SNR_start = 0;
SNR_step = 5;
SNR_end = 40;
SNR = SNR_start:SNR_step:SNR_end;% define SNR
SNR_lin = 10.^(SNR/10);% multipath linear SNR
C = [-3+3j -3+1j -3-3j -3-1j -1+3j -1+1j -1-3j -1-1j 1+3j 1+1j 1-3j 1-1j 3+3j 3+1j 3-3j 3-1j];% Graycode rule
Es=sum(abs(C).^2)/M;% average power of symbol
Eb=Es/symbol_bits;% average power of bit
N0=Eb./SNR_lin;% One-sided power spectral density
sigma=sqrt(N0/2);% standard deviation of AWGN
gamma_s=(Es.^2)./(sigma.^2);% compute for MMSE equalizer

%---MIMO parameter---
ntx = 2;% number of transmitters
nrx = [2 3 5 10];% number of receivers

%===Simulation===
Nav = 1e5;% number of average loops
%---define the size of variables---
C_Dk_zfe = zeros(1,ntx);
C_Dk_mmse = zeros(1,ntx);
rx_zfe = zeros(1,ntx);
rx_mmse = zeros(1,ntx);
ber_sim_16qam_mimo_zfe = zeros();
ber_sim_16qam_mimo_mmse = zeros();
for r=1:length(nrx)
    
    for m=1:length(SNR)
        
        %---Reset BER---
        total_errors_zfe = 0;
        total_errors_mmse = 0;
        
        
        for n=1:Nav
            %===Tx===
            %---16QAM modulation---
            tx = randi([0 M-1],1,ntx);% generate 16QAM symbols
            Dn_16qam_mod = C(tx+1).';% 16QAM modulation
            
            %===Channel===
            %---Generate MIMO channel matrix---
            H = (randn(nrx(r),ntx)+1j*randn(nrx(r),ntx))/sqrt(2);% Rayleigh Channel with complex Gaussian coefficients
            %---Compute C_zf---
            if ntx==nrx(r)
                C_zfe = inv(H);
            elseif ntx<nrx(r)
                C_zfe = pinv(H);
            end
            H_p = conj(H.');% Hermittian transpose operator for MMSE
            %---AWGN noise---
            awgn_noise = sigma(m)*(randn(nrx(r),1)+1i*randn(nrx(r),1));% generate AWGN noise
            signal_mimo = H*Dn_16qam_mod+awgn_noise*sqrt(ntx);% received signal
            
            %===Rx===
            %---Equalizer---
            C_mmse = inv(H_p*H+(1/gamma_s(m)))*H_p;% compute MMSE equalizer coefficient
            Dk_est_zfe = C_zfe*signal_mimo;% ZFE detection
            Dk_est_mmse=C_mmse*signal_mimo;% MMSE detection
            %---16QAM Demodulation---
            for l = 1:ntx
                C_Dk_zfe = sqrt(abs(C-Dk_est_zfe(l)).^2);% compute the Euclidean distance of every symbols to constellation points
                C_Dk_mmse = sqrt(abs(C-Dk_est_mmse(l)).^2);
                [b,index_zf] = min(C_Dk_zfe);% map  points to constellation points with shortest distance
                [c,index_mmse] = min(C_Dk_mmse);
                rx_zfe(l) = index_zf-1;% the original signal were 0-15
                rx_mmse(l) = index_mmse-1;
            end
            %===Compute BER===
            bit_errors_zf = biterr(tx,rx_zfe);% compare original signal and detected signal
            bit_errors_mmse = biterr(tx,rx_mmse);
            total_errors_zfe = total_errors_zfe+bit_errors_zf;% sum errors in average loops
            total_errors_mmse = total_errors_mmse+bit_errors_mmse;
        end
        ber_sim_16qam_mimo_zfe(m,r) = total_errors_zfe/symbol_bits/ntx/Nav;% compute BER
        ber_sim_16qam_mimo_mmse(m,r) = total_errors_mmse/symbol_bits/ntx/Nav;
        toc
    end
end

%===Plot results===
color_index=['bx-','bo-';'rx-' 'ro-';'gx-' 'go-';'kx-' 'ko-'];% colors used in graph
legend_index=cell(length(nrx)*2,1);% initializ legend set
for p=1:length(nrx)
    semilogy(SNR,ber_sim_16qam_mimo_zfe(:,p),color_index(p,1:3),SNR,ber_sim_16qam_mimo_mmse(:,p),color_index(p,4:6));% plot curves for a defined situation
    legend_index{p*2-1}=strcat(num2str(ntx),'x',num2str(nrx(p)),' ZF');% legends for ZF curves
    legend_index{p*2}=strcat(num2str(ntx),'x',num2str(nrx(p)),' MMSE');% legends for MMSE curves
    hold on;
end
legend(legend_index);
title('MIMO BER')
grid on
ylim([1e-6 1])
xlim([SNR_start SNR_end])
hold off;

%===Save results===
save_str=['MIMO_Sim_',num2str(Nav),'Nav'];
save(save_str);
toc