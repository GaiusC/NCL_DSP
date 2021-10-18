clear, clc, close all
tic
%===Simulation parameter definitions===
M=16;% M for 16 QAM
symbol_bits=log2(M);% 4 bits for a symbol in 16 QAM
nsymbol=2048;% Define 2048 symbols
ncarrier=2048;% Define 128 carriers
cp_length=16;% Length of Cyclic Prefix
SNR=0:1:42;% multipath SNR(dB)
SNR_lin=10.^(SNR/10);% multipath linear SNR
C=[-3+3j -3+1j -3-3j -3-1j -1+3j -1+1j -1-3j -1-1j 1+3j 1+1j 1-3j 1-1j 3+3j 3+1j 3-3j 3-1j];% Graycode rule
Es=sum(abs(C).^2)/M;% average power of symbol
Eb=Es/symbol_bits;
N0=Eb./SNR_lin;
sigma=sqrt(N0/2);
%---Channel Energy Normalization---
P=[0 -0.9 -4.9 -8 -7.8 -23.9];
h_l=sqrt(10.^(P/10));
U=sqrt(sum(abs(h_l).^2));
hl=h_l./U;
H=[hl(1),zeros(1,1),hl(2),zeros(1,6),hl(3),zeros(1,4),hl(4),zeros(1,11),hl(5),zeros(1,15),hl(6)];
Hk=fft(H,2048);

%===Semianalyical BER===
%---BPSK in multipath---
BER_BPSK_MP_symbol=zeros(1,ncarrier);
BER_BPSK_MP=zeros(1,length(SNR));
for mm=1:length(SNR)
    BER_BPSK_MP_symbol=qfunc(sqrt(2*SNR_lin(mm)*(abs(Hk).^2)));
    BER_BPSK_MP(mm)=sum(BER_BPSK_MP_symbol)/ncarrier;
end
%---BPSK AWGN---
BER_BPSK_AWGN=qfunc(sqrt(2*SNR_lin));
%---16QAM in multipath---
BER_16QAM_MP_symbol=zeros(1,ncarrier);
BER_16QAM_MP=zeros(1,length(SNR));
for kk=1:length(SNR)
    BER_16QAM_MP_symbol=qfunc(sqrt(3*symbol_bits/(M-1)*SNR_lin(kk)*abs(Hk).^2));
    BER_16QAM_MP(kk)=2*(1-1/sqrt(M))*sum(BER_16QAM_MP_symbol)/ncarrier;
end
Pm_MP=1-(1-BER_16QAM_MP).^2;
Pb_MP=Pm_MP./symbol_bits;
%---16QAM AWGN---
BER_16QAM_AWGN=zeros(1,length(SNR));
for oo=1:length(SNR)
    P_sqrtM = 2*(1-1./sqrt(M))*qfunc(sqrt(3*symbol_bits/(M-1)*SNR_lin(oo)));
    P_M = 1-(1-P_sqrtM).^2;
    BER_16QAM_AWGN(oo) = P_M/symbol_bits;
end

%===Simulation===
%---define the size of variables---
ber_sim_16qam_awgn=zeros(1,length(SNR));
ber_sim_16qam_mp=zeros(1,length(SNR));
zf=[];
for m=1:length(SNR)
    total_errors_awgn=0;
    total_errors_mp=0;
    %===Tx===
    %---16QAM modulation---
    tx=randi([0 M-1],1,nsymbol);
    Dn_16qam_mod=C(tx+1);
    %---ifft---
    Dn_time=ifft(Dn_16qam_mod)*sqrt(nsymbol);
    %---Insert Cyclic Prefix---
    xn=[Dn_time(1,nsymbol-cp_length+1:end),Dn_time];
    
    %===Channel===
    %---AWGN Config---
    [yn,zf]=filter(H,1,xn,zf);
    awgn_noise=sigma(m)*(randn(1,nsymbol+cp_length)+1i*randn(1,nsymbol+cp_length));
    %---Add AWGN---
    signal_awgn=xn+awgn_noise;
    signal_mp=yn+awgn_noise;
    
    %===Rx===
    %---Remove Cyclic Prefix---
    rx_cp_awgn=signal_awgn(:,cp_length+1:end);
    rx_cp_mp=signal_mp(:,cp_length+1:end);
    %---fft---
    rx_signal_freq_awgn=fft(rx_cp_awgn)/sqrt(nsymbol);
    rx_signal_freq_mp=fft(rx_cp_mp)/sqrt(nsymbol);
    %---16QAM Demodulation---
    C_Dk_awgn=zeros(1,M);
    C_Dk_mp=zeros(1,M);
    rx_awgn=zeros(1,nsymbol);
    rx_mp=zeros(1,nsymbol);
    for l =1:nsymbol
        C_Dk_awgn = sqrt(abs(C-rx_signal_freq_awgn(l)).^2);
        C_Dk_mp = sqrt(abs(C-rx_signal_freq_mp(l)).^2);
        [a,index_awgn] = min(C_Dk_awgn);
        [b,index_mp] = min(C_Dk_mp);
        rx_awgn(l) = index_awgn-1;
        rx_mp(l) = index_mp-1;
    end
    
    %===Compute BER===
    bit_errors_awgn=biterr(rx_awgn,tx);
    bit_errors_mp=biterr(rx_mp,tx);
    total_errors_awgn=total_errors_awgn+bit_errors_awgn;
    total_errors_mp=total_errors_mp+bit_errors_mp;
    ber_sim_16qam_awgn(m)=total_errors_awgn/nsymbol/symbol_bits;
    ber_sim_16qam_mp(m)=total_errors_mp/nsymbol/symbol_bits;
    toc
end

%===plot results===
figure()
semilogy(SNR,BER_BPSK_AWGN,'c-',SNR,BER_BPSK_MP,'r-',SNR,BER_16QAM_AWGN,'b-',SNR,ber_sim_16qam_awgn,'gs',SNR,Pb_MP,'m-',SNR,ber_sim_16qam_mp,'x-k');
legend('BPSK AWGN','BPSK MP','16 QAM AWGN(Theory)','16 QAM MP(Semi.)','16 QAM MF MP (Sim.)');
ylim([1e-5 1]);
xlim([0 30]);
grid on;
