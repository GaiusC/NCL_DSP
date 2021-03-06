%Author: J.Chen
%Matlab Version: R2021a
clear, clc, close all
tic
%===Reset random number===
RN=sum(100*clock);
RS=RandStream('mt19937ar','seed',RN);
RandStream.setGlobalStream(RS);
%===Simulation parameter definitions===
M=16;% M for 16 QAM
symbol_bits=log2(M);% 4 bits for a symbol in 16 QAM
ncarriers=2048;% Define 2048 symbols
cp_length=64;% Length of Cyclic Prefix
SNR_mp=0:1:30;% multipath SNR(dB)
SNR_mp_lin=10.^(SNR_mp/10);% multipath linear SNR
C=[-3+3j -3+1j -3-3j -3-1j -1+3j -1+1j -1-3j -1-1j 1+3j 1+1j 1-3j 1-1j 3+3j 3+1j 3-3j 3-1j];% Graycode rule
Es=sum(abs(C).^2)/M;% average power of symbol
Eb=Es/symbol_bits;% average power of bit
N0=Eb./SNR_mp_lin;% One-sided power spectral density
sigma=sqrt(N0/2);% standard deviation of AWGN
gama_s=(Es.^2)./(sigma.^2);% compute for MMSE equalizer

%===Channel Energy Normalization===
P=[0 -0.9 -4.9 -8 -7.8 -23.9];% define channel for scenario
h_l=sqrt(10.^(P/10));% channel magnitudes
U=sqrt(sum(abs(h_l).^2));% normalize coefficient
hl=h_l./U;% normalize coefficient
H=[hl(1),zeros(1,1),hl(2),zeros(1,6),hl(3),zeros(1,4),hl(4),zeros(1,11),hl(5),zeros(1,15),hl(6)];% Zero-padding
Hk=fft(H,ncarriers);% channel frequency response (CFR)
C_zf=1./Hk;% ZF equalizer coefficient

%===Theoretical BER in Multipath===
%---16QAM in multipath---
BER_16QAM_MP=zeros(1,length(SNR_mp));
for kk=1:length(SNR_mp)
    BER_16QAM_MP_symbol=qfunc(sqrt(3*symbol_bits/(M-1)*SNR_mp_lin(kk)*abs(Hk).^2));
    BER_16QAM_MP(kk)=2*(1-1/sqrt(M))*sum(BER_16QAM_MP_symbol)/ncarriers;
end
Pm_MP=1-(1-BER_16QAM_MP).^2;
Pb_MP=Pm_MP./symbol_bits;

%---BPSK in multipath---
BER_BPSK_MP_symbol=zeros(1,ncarriers);
BER_BPSK_MP=zeros(1,length(SNR_mp));
for mm=1:length(SNR_mp)
    BER_BPSK_MP_symbol=qfunc(sqrt(2*SNR_mp_lin(mm)*(abs(Hk).^2)));
    BER_BPSK_MP(mm)=sum(BER_BPSK_MP_symbol)/ncarriers;
end

%---BPSK AWGN---
BER_BPSK_AWGN=qfunc(sqrt(2*SNR_mp_lin));

%---16QAM AWGN---
BER_16QAM_AWGN=zeros(1,length(SNR_mp));
for oo=1:length(SNR_mp)
    P_sqrtM = 2*(1-1./sqrt(M))*qfunc(sqrt(3*symbol_bits/(M-1)*SNR_mp_lin(oo)));
    P_M = 1-(1-P_sqrtM).^2;
    BER_16QAM_AWGN(oo) = P_M/symbol_bits;
end
toc

%===Simulation===
%---define the size of variables---
C_Dk_mf=zeros(1,M);
C_Dk_awgn=zeros(1,M);
C_Dk_zf=zeros(1,M);
C_Dk_mmse=zeros(1,M);
rx_mf=zeros(1,ncarriers);
rx_awgn=zeros(1,ncarriers);
rx_zf=zeros(1,M);
rx_mmse=zeros(1,M);
ber_sim_16qam_mf=zeros(1,length(SNR_mp));
ber_sim_16qam_mp_awgn=zeros(1,length(SNR_mp));
ber_sim_16qam_mp_zf=zeros(1,length(SNR_mp));
ber_sim_16qam_mp_mmse=zeros(1,length(SNR_mp));
zf=[];
Nav=500;% number of average loops
%---loop---
for m=1:length(SNR_mp)
    %---Reset---
    total_errors_mf=0;
    total_errors_awgn=0;
    total_errors_zf=0;
    total_errors_mmse=0;
    for n=1:Nav
        
        %===Tx===
        %---16QAM modulation---
        tx=randi([0 M-1],1,ncarriers);% generate 16QAM symbols
        Dn_16qam_mod=C(tx+1);% 16QAM modulation
        %---ifft---
        Dn=ifft(Dn_16qam_mod)*sqrt(ncarriers);% convert signal to time-series and energy normalization
        %---Insert Cyclic Prefix---
        xn=[Dn(1,ncarriers-cp_length+1:end),Dn];% use the last 64 data of Dn to be its own CP
        
        %===Channel===
        %---AWGN Config---
        [yn,zf]=filter(H,1,xn,zf);% convolution the signal with channel
        awgn_noise=sigma(m)*(randn(1,ncarriers+cp_length)+1i*randn(1,ncarriers+cp_length));% generate AWGN noise
        %---Add AWGN---
        signal_mf=xn+awgn_noise;% received for simulation in sigle path
        signal_awgn=yn+awgn_noise;% received for simulation with no equalizer
        signal_known=yn+awgn_noise;% received for simulation with 2 equalizers
        
        %===Rx===
        %---Remove Cyclic Prefix---
        rx_cp_mf=signal_mf(:,cp_length+1:end);
        rx_cp_awgn=signal_awgn(:,cp_length+1:end);
        rx_cp_known=signal_known(:,cp_length+1:end);
        %---fft---
        rx_signal_freq_mf=fft(rx_cp_mf)/sqrt(ncarriers);% convert signal to frequency-series and energy normalization
        rx_signal_freq_awgn=fft(rx_cp_awgn)/sqrt(ncarriers);
        rx_signal_freq_known=fft(rx_cp_known)/sqrt(ncarriers);
        %---One-tap equalizer---
        C_mmes=conj(Hk)./((abs(Hk).^2)+(1/gama_s(m)));% MMSE equalizer coefficient
        Dk_est_zf=C_zf.*rx_signal_freq_known;% ZFE detection
        Dk_est_mmse=C_mmes.*rx_signal_freq_known;% MMSE detection
        %---16QAM Demodulation---
        for l =1:ncarriers
            C_Dk_mf = sqrt(abs(C-rx_signal_freq_mf(l)).^2);% compute the Euclidean distance of every symbols to constellation points
            C_Dk_awgn = sqrt(abs(C-rx_signal_freq_awgn(l)).^2);
            C_Dk_zf = sqrt(abs(C-Dk_est_zf(l)).^2);
            C_Dk_mmse = sqrt(abs(C-Dk_est_mmse(l)).^2);
            [a,index_mf] = min(C_Dk_mf);% map  points to constellation points with shortest distance
            [b,index_awgn] = min(C_Dk_awgn);
            [c,index_zf] = min(C_Dk_zf);
            [d,index_mmse] = min(C_Dk_mmse);
            rx_mf(l) = index_mf-1;% the original signal were 0-15
            rx_awgn(l) = index_awgn-1;
            rx_zf(l) = index_zf-1;
            rx_mmse(l) = index_mmse-1;
        end
        %===Compute BER===
        bit_errors_mf=biterr(tx,rx_mf);% compute the Euclidean distance of every symbols to constellation points
        bit_errors_awgn=biterr(tx,rx_awgn);
        bit_errors_zf=biterr(tx,rx_zf);
        bit_errors_mmse=biterr(tx,rx_mmse);
        total_errors_mf=total_errors_mf+bit_errors_mf;% sum errors in average loops
        total_errors_awgn=total_errors_awgn+bit_errors_awgn;
        total_errors_zf=total_errors_zf+bit_errors_zf;
        total_errors_mmse=total_errors_mmse+bit_errors_mmse;
    end
    ber_sim_16qam_mf(m)=total_errors_mf/ncarriers/symbol_bits/Nav;% compute BER
    ber_sim_16qam_mp_awgn(m)=total_errors_awgn/ncarriers/symbol_bits/Nav;
    ber_sim_16qam_mp_zf(m)=total_errors_zf/ncarriers/symbol_bits/Nav;
    ber_sim_16qam_mp_mmse(m)=total_errors_mmse/ncarriers/symbol_bits/Nav;
    
    toc
end

%plot results
figure()
semilogy(SNR_mp,BER_BPSK_AWGN,'c-',SNR_mp,BER_BPSK_MP,'r-',SNR_mp,BER_16QAM_AWGN,'b-',SNR_mp,ber_sim_16qam_mf,'gs',SNR_mp,Pb_MP,'m-',SNR_mp,ber_sim_16qam_mp_awgn,'x-k',SNR_mp,ber_sim_16qam_mp_zf,'o-g',SNR_mp,ber_sim_16qam_mp_mmse,'*-y');
legend('BPSK AWGN','BPSK MP','16 QAM AWGN(Theory)','16 QAM MF(Sim.)','16 QAM MP(Semi.)','16 QAM MF MP (Sim.)','16 QAM ZF MP(Sim.)','16 QAM MMSE MP(Sim.)');
ylim([1e-5 1]);
xlim([0 30]);
grid on;

%===Save results===
save MP_KnownCH;