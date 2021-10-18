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
SNR_step = 1;
SNR_end = 40;
SNR = SNR_start:SNR_step:SNR_end;
SNR_lin = 10.^(SNR/10);% multipath linear SNR
C = [-3+3j -3+1j -3-3j -3-1j -1+3j -1+1j -1-3j -1-1j 1+3j 1+1j 1-3j 1-1j 3+3j 3+1j 3-3j 3-1j];% Graycode rule
Es=sum(abs(C).^2)/M;% average power of symbol
Eb=Es/symbol_bits;
N0=Eb./SNR_lin;
sigma=sqrt(N0/2);
gamma_s=(Es.^2)./(sigma.^2);

%---MIMO parameter---
ntx = 2;
nrx = 2;

%===Simulation===
Nav = 1e4;
%---define the size of variables---
C_Dk_zfe = zeros(1,ntx);
C_Dk_mmse = zeros(1,ntx);
rx_zfe = zeros(1,ntx);
rx_mmse = zeros(1,ntx);
ber_sim_16qam_mimo_zfe = zeros();
ber_sim_16qam_mimo_mmse = zeros();

for m=1:length(SNR)
    
    %---Reset---
    total_errors_zfe = 0;
    total_errors_mmse = 0;
    
    
    for n=1:Nav
        %===Tx===
        %---16QAM modulation---
        tx = randi([0 M-1],1,ntx);
        Dn_16qam_mod = C(tx+1).';
        
        %===Channel===
        %---Generate MIMO channel matrix---
        H = (randn(nrx,ntx)+1j*randn(nrx,ntx))/sqrt(2);
        %---Compute C_zf---
        if ntx==nrx
            C_zfe = inv(H);
        elseif ntx<nrx
            C_zfe = pinv(H);
        end
        H_p = conj(H.');
        awgn_noise = sigma(m)*(randn(nrx,1)+1i*randn(nrx,1));
        signal_mimo = H*Dn_16qam_mod+awgn_noise*sqrt(ntx);
        
        %===Rx===
        %---Equalizer---
        C_mmse = inv(H_p*H+(1/gamma_s(m)))*H_p;%
        Dk_est_zfe = C_zfe*signal_mimo;
        Dk_est_mmse=C_mmse*signal_mimo;
        %---16QAM Demodulation---
        for l = 1:ntx
            C_Dk_zfe = sqrt(abs(C-Dk_est_zfe(l)).^2);
            C_Dk_mmse = sqrt(abs(C-Dk_est_mmse(l)).^2);
            [b,index_zf] = min(C_Dk_zfe);
            [c,index_mmse] = min(C_Dk_mmse);
            rx_zfe(l) = index_zf-1;
            rx_mmse(l) = index_mmse-1;
        end
        
        %===Compute BER===
        bit_errors_zf = biterr(tx,rx_zfe);
        bit_errors_mmse = biterr(tx,rx_mmse);
        total_errors_zfe = total_errors_zfe+bit_errors_zf;
        total_errors_mmse = total_errors_mmse+bit_errors_mmse;
    end
    ber_sim_16qam_mimo_zfe(m) = total_errors_zfe/symbol_bits/ntx/Nav;
    ber_sim_16qam_mimo_mmse(m) = total_errors_mmse/symbol_bits/ntx/Nav;
    
    %===Print===
    str=sprintf('SNR=%f. BER_zfe=%f, BER_mmse=%f/n',SNR(m),ber_sim_16qam_mimo_zfe(m),ber_sim_16qam_mimo_mmse(m));
    disp(str);
    toc
end

%===Plot results===
semilogy(SNR,ber_sim_16qam_mimo_zfe,'r',SNR,ber_sim_16qam_mimo_mmse,'b')
title('MIMO BER')
legend('ZFE','MMSE')
grid on
ylim([1e-5 1])
xlim([SNR_start SNR_end])

true(ber_sim_16qam_mimo_mmse(1)<ber_sim_16qam_mimo_zfe(1))

toc