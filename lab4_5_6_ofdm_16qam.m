clear, clc, close all
%===Simulation parameter definitions===
M=16; % Constellation size
nsymbol=2048;% Define 2048 symbols
symbol_bits=log2(M);% 4 bits per symbol
cp_length=64;% Length of Cyclic Prefix
C=[-3+3j -3+1j -3-3j -3-1j -1+3j -1+1j -1-3j -1-1j 1+3j 1+1j 1-3j 1-1j 3+3j 3+1j 3-3j 3-1j];% Graycode rule
spow=mean(abs(C).^2);% average power of symbols
SNR=0:1:30;% Range of SNR
SNR_lin=10.^(SNR/10);% Convert dB to linear value
sigma=sqrt(0.5.*spow./symbol_bits./SNR_lin);

%===theoretical BER===
for n=1:length(SNR)
    P_sqrtM = 2*(1-1./sqrt(M))*qfunc(sqrt(3*symbol_bits/(M-1)*SNR_lin(n)));
    P_M = 1-(1-P_sqrtM).^2;
    P_b(n) = P_M/symbol_bits;
end

%===simulation===
%---define size---
rx=zeros(1,nsymbol);
BER=zeros(length(SNR),nsymbol);
for m=1:length(SNR)
    
    total_errors=0;
    %==========Tx==========
    %---16QAM modulation---
    tx=randi([0 M-1],1,nsymbol);% Generate random order in signal
    Dn_16qam_mod=C(tx+1);% Map the order to generate 16QAM signal
    %---ifft---
    tx_signal_time=ifft(Dn_16qam_mod)*sqrt(nsymbol);
    %---Insert Cyclic Prefix---
    dn_cp=[tx_signal_time(:,1:cp_length),tx_signal_time];
    
    %=======Channel========
    awgn_noise=sigma(m)*(randn(1,nsymbol+cp_length)+1i*randn(1,nsymbol+cp_length));
    signal=dn_cp+awgn_noise;
    
    %==========Rx==========
    %---Remove Cyclic Prefix---
    rx_cp=signal(1,cp_length+1:end);
    %---fft---
    rx_signal_freq=fft(rx_cp)/sqrt(nsymbol);
    %---16QAM Demodulation---
    for l =1:nsymbol
        C_Dk = sqrt(abs(C-rx_signal_freq(l)).^2);
        [a,index] = min(C_Dk);
        rx(l) = index-1;
    end
    
    %======Compute BER======
    new_errors=biterr(rx,tx);
    total_errors=total_errors+new_errors;
    BER(m)=total_errors/nsymbol/symbol_bits;
end

%======Plot Results======
scatterplot(Dn_16qam_mod);
scatterplot(rx_signal_freq);
figure()
semilogy(SNR,BER,"ro",SNR,P_b,"b-");
ylim([1e-5 1]);
xlim([0 15]);
title("BER of 16-QAM in OFDM")
xlabel("Eb/N0 (dB)");
ylabel("BER");
legend("BER", "Theoretical BER");
grid on;

%======Save Results======
save lab4_5_6_ofdm_16qam_results
