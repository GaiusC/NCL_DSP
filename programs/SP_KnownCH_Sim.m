clear, clc, close all
tic
%===Simulation===
%---BPSK---
nsymbol_BPSK=1e6;% Define 100000 symbols
SNR_BPSK=-5:5:45;% multipath SNR(dB)
SNR_lin_BPSK=10.^(SNR_BPSK/10);% multipath linear SNR
C_BPSK=[1 -1];% BPSK code
BER_sim_BPSK=zeros(length(SNR_BPSK),1);
dn_r_BPSK=zeros(1,nsymbol_BPSK);
bn_r_BPSK=zeros(1,nsymbol_BPSK);

for m=1:length(SNR_BPSK)
    
    total_errors_BPSK=0;
    total_errors_16QAM=0;
    
    %Tx
    tx_BPSK=randi([0 1],1,nsymbol_BPSK);% generate BPSK signals
    Dn_BPSK=C_BPSK(tx_BPSK+1);% BPSK modulation
    
    %Channel
    hn=1/sqrt(2)*(randn(1,nsymbol_BPSK)+1i*randn(1,nsymbol_BPSK));% generate Rayleigh Channel
    sigma_BPSK=sqrt(0.5/SNR_lin_BPSK(m));% compute standard deviation for AWGN
    wn=sigma_BPSK*(randn(1,nsymbol_BPSK)+1i*randn(1,nsymbol_BPSK));% AWGN noise
    rn=hn.*Dn_BPSK+wn;% received signal
    %Rx
    for n=1:nsymbol_BPSK
        dn_r_BPSK(n)=real(rn(n)./hn(n));% BPSK demodulation
        if dn_r_BPSK(n)<=0
            bn_r_BPSK(n)=1;
        else
            bn_r_BPSK(n)=0;
        end
    end
    %BER
    new_errors=biterr(bn_r_BPSK,tx_BPSK);% compare original signal and detected signal
    total_errors_BPSK=total_errors_BPSK+new_errors;% sum errors in average loops
    BER_sim_BPSK(m)=total_errors_BPSK/nsymbol_BPSK;% compute BER
end
toc
%---16QAM---
M=16;% M for 16 QAM
nsymbol_16QAM=2048;% Define 2048 symbols
symbol_bits=log2(M);% 4 bits per symbol
cp_length=64;% Length of Cyclic Prefix
C_16QAM=[-3+3j -3+1j -3-3j -3-1j -1+3j -1+1j -1-3j -1-1j 1+3j 1+1j 1-3j 1-1j 3+3j 3+1j 3-3j 3-1j];% Graycode rule
spow=mean(abs(C_16QAM).^2);% average power of symbols
SNR_16QAM=0:1:30;% Range of SNR
SNR_lin_16QAM=10.^(SNR_16QAM/10);% Convert dB to linear value
sigma_16QAM=sqrt(0.5.*spow./symbol_bits./SNR_lin_16QAM);% standard deviation of AWGN
BER_sim_16QAM=zeros(length(SNR_16QAM),nsymbol_16QAM);

Nav=500;% number of average loops
for m=1:length(SNR_16QAM)
    total_errors_16QAM=0;% reset BER
    for n=1:Nav
        
        %==========Tx==========
        %---16QAM modulation---
        tx_16QAM=randi([0 M-1],1,nsymbol_16QAM);% Generate random order in signal
        Dn_16qam_mod=C_16QAM(tx_16QAM+1);% Map the order to generate 16QAM signal
        %---ifft---
        tx_signal_time=ifft(Dn_16qam_mod)*sqrt(nsymbol_16QAM);% convert signal to time-series
        %---Insert Cyclic Prefix---
        dn_cp=[tx_signal_time(:,1:cp_length),tx_signal_time];% use the last n data to be the CP
        
        %=======Channel========
        awgn_noise_025=0.25*(randn(1,nsymbol_16QAM)+1i*randn(1,nsymbol_16QAM));% AWGN noise for standard deviation as 0.25
        awgn_noise=sigma_16QAM(m)*(randn(1,nsymbol_16QAM+cp_length)+1i*randn(1,nsymbol_16QAM+cp_length));% generate AWGN noise
        signal025=Dn_16qam_mod+awgn_noise_025;% received signal with low standard deviation
        signal=dn_cp+awgn_noise;% received signal
        
        %==========Rx==========
        %---Remove Cyclic Prefix---
        rx_cp=signal(1,cp_length+1:end);
        %---fft---
        rx_signal_freq=fft(rx_cp)/sqrt(nsymbol_16QAM);% convert signal to frequency-series
        %---16QAM Demodulation---
        for l =1:nsymbol_16QAM
            C_Dk = sqrt(abs(C_16QAM-rx_signal_freq(l)).^2);% compute the Euclidean distance of every symbols to constellation points
            [a,index] = min(C_Dk);% map  points to constellation points with shortest distance
            rx(l) = index-1;% the original signal were 0-15
        end
        
        %======Compute BER======
        new_errors=biterr(rx,tx_16QAM);% compare original signal and detected signal
        total_errors_16QAM=total_errors_16QAM+new_errors;% sum errors in average loops
    end
    BER_sim_16QAM(m)=total_errors_16QAM/nsymbol_16QAM/symbol_bits/Nav;% compute BER
    toc
end

%===Theoretical BER===
%---BPSK in AWGN---
gamma_dB=-5:1:10;
SNR_lin_BPSK=10.^(gamma_dB/10);
ber_AWGN=qfunc(sqrt(2*SNR_lin_BPSK));

%---BPSK in Rayleigh---
Gamma_dB=-5:5:50;
Gamma_lin=10.^(Gamma_dB/10);
ber_Rayleigh=(1-sqrt(Gamma_lin./(Gamma_lin+1)))/2;

%===Semi-analytical===
%---BPSK in Rayleigh---
N=1e6;
gamma_dB_2=-5:5:50;
SNR_lin_BPSK=10.^(gamma_dB_2/10);
L=length(gamma_dB_2);
ber_Semi_Rayleigh=zeros(1,L);
for i=1:L
    h=1/sqrt(2)*(randn(1,N)+1j*randn(1,N));
    Pe=qfunc(sqrt(2*SNR_lin_BPSK(i)*abs(h).^2));
    ber_Semi_Rayleigh(i)=sum(Pe)/N;
    toc
end

%---16QAM in AWGN---
for n=1:length(SNR_16QAM)
    P_sqrtM = 2*(1-1./sqrt(M))*qfunc(sqrt(3*symbol_bits/(M-1)*SNR_lin_16QAM(n)));
    P_M = 1-(1-P_sqrtM).^2;
    P_b(n) = P_M/symbol_bits;
end

%===Plot results===
%---BPSK---
semilogy(gamma_dB,ber_AWGN,'r-',Gamma_dB,ber_Rayleigh,'b--',gamma_dB_2,ber_Semi_Rayleigh,'o--',SNR_BPSK,BER_sim_BPSK,'go--');
xlabel("Average SNR ");
ylabel("Pe");
ylim([1e-6 1]);
title('BER of BPSK in Different Channel');
grid on;
legend("BPSK in AWGN","BPSK in Rayleigh","Semi-analytial","Full Simulation");
%---16QAM---
figure();
plot(signal025,'ro');
hold on;
plot(Dn_16qam_mod,'ko','MarkerFaceColor','k');
hold off;
title('16QAM constellation in low standard deviation(0.25)');
grid on;
figure()
semilogy(SNR_16QAM,BER_sim_16QAM,"ro",SNR_16QAM,P_b,"b-");
ylim([1e-5 1]);
xlim([0 15]);
title("BER of 16-QAM in OFDM");
xlabel("Eb/N0 (dB)");
ylabel("BER");
legend("BER", "Theoretical BER");
grid on;
toc
%save results
save SP_KnownCH;
