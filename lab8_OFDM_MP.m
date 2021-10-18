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
Eb=Es/symbol_bits;
N0=Eb./SNR_mp_lin;
sigma=sqrt(N0/2);
gama_s=(Es.^2)./(sigma.^2);

%===Channel Energy Normalization===
P=[0 -0.9 -4.9 -8 -7.8 -23.9];
h_l=sqrt(10.^(P/10));
U=sqrt(sum(abs(h_l).^2));
hl=h_l./U;
H=[hl(1),zeros(1,1),hl(2),zeros(1,6),hl(3),zeros(1,4),hl(4),zeros(1,11),hl(5),zeros(1,15),hl(6)];
Hk=fft(H,ncarriers);
C_zf=1./Hk;

%===Theoretical BER in Multipath===
%---16QAM in multipath---
BER_16QAM_MP_symbol=zeros(1,ncarriers);
BER_16QAM_MP=zeros(1,length(SNR_mp));
for kk=1:length(SNR_mp)
    BER_16QAM_MP_symbol=qfunc(sqrt(3*symbol_bits/(M-1)*SNR_mp_lin(kk)*abs(Hk).^2));
    BER_16QAM_MP(kk)=2*(1-1/sqrt(M))*sum(BER_16QAM_MP_symbol)/ncarriers;
end
Pm_MP=1-(1-BER_16QAM_MP).^2;
Pb_MP=Pm_MP./symbol_bits;

%===Simulation===
%---define the size of variables---
C_Dk_awgn=zeros(1,M);
C_Dk_zf=zeros(1,M);
C_Dk_mmse=zeros(1,M);
rx_awgn=zeros(1,ncarriers);
rx_zf=zeros(1,M);
rx_mmse=zeros(1,M);
ber_sim_16qam_mp_awgn=zeros(1,length(SNR_mp));
ber_sim_16qam_mp_zf=zeros(1,length(SNR_mp));
ber_sim_16qam_mp_mmse=zeros(1,length(SNR_mp));
zf=[];

Nav=100;
%---loop---
for m=1:length(SNR_mp)
    %---Reset---
    total_errors_awgn=0;
    total_errors_zf=0;
    total_errors_mmse=0;
    for n=1:Nav
        
        %===Tx===
        %---16QAM modulation---
        tx=randi([0 M-1],1,ncarriers);
        Dn_16qam_mod=C(tx+1);
        %---ifft---
        Dn=ifft(Dn_16qam_mod)*sqrt(ncarriers);% '*sqrt(nsymbol)' for energy normalization
        %---Insert Cyclic Prefix---
        xn=[Dn(1,ncarriers-cp_length+1:end),Dn];% use the last 64 data of Dn to be its own CP
        
        %===Channel===
        %---AWGN Config---
        [yn,zf]=filter(H,1,xn,zf);
        awgn_noise=sigma(m)*(randn(1,ncarriers+cp_length)+1i*randn(1,ncarriers+cp_length));
        %---Add AWGN---
        signal_awgn=yn+awgn_noise;
        signal_known=yn+awgn_noise;
        
        %===Rx===
        %---Remove Cyclic Prefix---
        rx_cp_awgn=signal_awgn(:,cp_length+1:end);
        rx_cp_known=signal_known(:,cp_length+1:end);
        %---fft---
        rx_signal_freq_awgn=fft(rx_cp_awgn)/sqrt(ncarriers);
        rx_signal_freq_known=fft(rx_cp_known)/sqrt(ncarriers);
        %---One-tap equalizer---
        %C_zf=conj(Hk)./(abs(Hk).^2);
        C_mmes=conj(Hk)./((abs(Hk).^2)+(1/gama_s(m)));
        Dk_est_zf=C_zf.*rx_signal_freq_known;
        Dk_est_mmse=C_mmes.*rx_signal_freq_known;
        %---16QAM Demodulation---
        for l =1:ncarriers
            C_Dk_awgn = sqrt(abs(C-rx_signal_freq_awgn(l)).^2);
            C_Dk_zf = sqrt(abs(C-Dk_est_zf(l)).^2);
            C_Dk_mmse = sqrt(abs(C-Dk_est_mmse(l)).^2);
            [a,index_awgn] = min(C_Dk_awgn);
            [b,index_zf] = min(C_Dk_zf);
            [c,index_mmse] = min(C_Dk_mmse);
            rx_awgn(l) = index_awgn-1;
            rx_zf(l) = index_zf-1;
            rx_mmse(l) = index_mmse-1;
        end
        %===Compute BER===
        bit_errors_awgn=biterr(tx,rx_awgn);
        bit_errors_zf=biterr(tx,rx_zf);
        bit_errors_mmse=biterr(tx,rx_mmse);
        total_errors_awgn=total_errors_awgn+bit_errors_awgn;
        total_errors_zf=total_errors_zf+bit_errors_zf;
        total_errors_mmse=total_errors_mmse+bit_errors_mmse;
    end
    
    ber_sim_16qam_mp_awgn(m)=total_errors_awgn/ncarriers/symbol_bits/Nav;
    ber_sim_16qam_mp_zf(m)=total_errors_zf/ncarriers/symbol_bits/Nav;
    ber_sim_16qam_mp_mmse(m)=total_errors_mmse/ncarriers/symbol_bits/Nav;
    
    toc
end

%plot results
figure()
semilogy(SNR_mp,Pb_MP,'b-',SNR_mp,ber_sim_16qam_mp_awgn,'x-k',SNR_mp,ber_sim_16qam_mp_zf,'o-g',SNR_mp,ber_sim_16qam_mp_mmse,'*-y');
legend('16 QAM MP(Semi.)','16 QAM MF MP (Sim.)','16 QAM ZF MP(Sim.)','16 QAM MMSE MP(Sim.)');
ylim([1e-5 1]);
xlim([0 30]);
grid on;
