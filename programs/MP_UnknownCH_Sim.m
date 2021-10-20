%Author: J.Chen
%Matlab Version: R2021a
clear, clc, close all
tic
%===Reset random number===
RN=sum(100*clock);
RS=RandStream('mt19937ar','seed',RN);
RandStream.setGlobalStream(RS);

%===Simulation parameter definitions===
M = 16;% M for 16 QAM
symbol_bits = log2(M);% 4 bits for a symbol in 16 QAM
ncarriers = 2048;% Define 2048 carriers
cp_length = [0 8 16 32 64];% Length of Cyclic Prefix
SNR_mp = 0:1:30;% multipath SNR(dB)
SNR_mp_lin = 10.^(SNR_mp/10);% multipath linear SNR
C = [-3+3j -3+1j -3-3j -3-1j -1+3j -1+1j -1-3j -1-1j 1+3j 1+1j 1-3j 1-1j 3+3j 3+1j 3-3j 3-1j];% Graycode rule
Es=sum(abs(C).^2)/M;% average power of symbol
Eb=Es/symbol_bits;% average power of bit
N0=Eb./SNR_mp_lin;% One-sided power spectral density
sigma=sqrt(N0/2);% standard deviation of AWGN
gamma_s=(Es.^2)./(sigma.^2);% compute for MMSE equalizer
pilot_content = 3+3i;% pilot content
pilot_position = 1:8:2048;% pilot position

%===Channel Energy Normalization===
P=[0 -0.9 -4.9 -8 -7.8 -23.9];% define channel for scenario
h_l=sqrt(10.^(P/10));% channel magnitudes
U=sqrt(sum(abs(h_l).^2));% normalize coefficient
hl=h_l./U;% normalize coefficient
H=[hl(1),zeros(1,1),hl(2),zeros(1,6),hl(3),zeros(1,4),hl(4),zeros(1,11),hl(5),zeros(1,15),hl(6)];% Zero-padding
Hk=fft(H,ncarriers);% channel frequency response (CFR)


%===Theoretical BER===
%---16QAM in multipath---
BER_16QAM_MP_symbol=zeros(1,ncarriers);
BER_16QAM_MP=zeros(1,length(SNR_mp));
for kk=1:length(SNR_mp)
    BER_16QAM_MP_symbol=qfunc(sqrt(3*symbol_bits/(M-1)*SNR_mp_lin(kk)*abs(Hk).^2));
    BER_16QAM_MP(kk)=2*(1-1/sqrt(M))*sum(BER_16QAM_MP_symbol)/ncarriers;
end
Pm_MP=1-(1-BER_16QAM_MP).^2;
Pb_MP=Pm_MP./symbol_bits;

toc
%===define the size of variables===
Nav = 500;% number of average loops
C_Dk_zf = zeros(1,M);
C_Dk_mmse = zeros(1,M);
Hk_est = zeros(1,ncarriers);
rx_zf = zeros(1,M);
rx_mmse = zeros(1,M);
ber_sim_16qam_mp_zf = zeros(1,length(SNR_mp));
ber_sim_16qam_mp_mmse = zeros(1,length(SNR_mp));
zf = [];
MSE_symbol = zeros(1,Nav);
MSE = zeros(length(SNR_mp),length(cp_length));

for m = 1:length(SNR_mp)
    
    for n = 1:length(cp_length)
        
        total_errors_zf = 0;
        total_errors_mmse = 0;
        
        for nav=1:Nav
            
            %===Tx===
            %---16QAM modulation---
            tx = randi([0 M-1],1,ncarriers);% generate 16QAM symbols
            Dn_16qam_mod = C(tx+1);% 16QAM modulation
            %---Insert pilot---
            Dn_pilot = Dn_16qam_mod;
            Dn_pilot(1,pilot_position) = pilot_content;% overwrite the original signal at the pilot position
            %---ifft---
            signal_tx = ifft(Dn_pilot)*sqrt(ncarriers);% convert signal to time-series and energy normalization
            %---cp---
            xn = [signal_tx(1,ncarriers-cp_length(n)+1:end),signal_tx];% use the last n data to be the CP
            
            %===Channel===
            [yn,zf] = filter(H,1,xn,zf);% convolution the signal with channel
            awgn_noise = sigma(m)*(randn(1,ncarriers+cp_length(n))+1i*randn(1,ncarriers+cp_length(n)));% generate AWGN noise
            signal_awgn = yn+awgn_noise;% received signal for simulation with equalizers
            
            %===Rx===
            %---Remove Cyclic Prefix---
            rx_cp = signal_awgn(:,cp_length(n)+1:end);
            %---fft---
            signal_rx = fft(rx_cp)/sqrt(ncarriers);% convert signal to frequency-series and energy normalization
            %---equalizer---
            Hp_est = signal_rx(pilot_position)/pilot_content;% estimate the channel coefficient using known pilot content and position
            k = 0;
            for pp = 1:255% there are 255 pilots
                for ii = 0:7% 7 subcarriers between 2 pilots
                    k = k+1;
                    Hk_est(k) = Hp_est(pp)+(ii/8)*(Hp_est(pp+1)-Hp_est(pp));% linear interpolation
                end
            end
            for iii = 8:15% estimate the last 7 index
                k = k+1;
                Hk_est(k) = Hp_est(pp)+(iii/8)*(Hp_est(pp+1)-Hp_est(pp));% linear interpolation
            end
            C_zf = conj(Hk_est)./(abs(Hk_est).^2);% ZF equalizer coefficient
            C_mmse = conj(Hk_est)./((abs(Hk_est).^2)+(1./gamma_s(m)));% MMSE equalizer coefficient
            Dk_est_zf = C_zf.*signal_rx;% ZFE detection
            Dk_est_zf(:,pilot_position) = [];% eliminate data on pilot position
            Dk_est_mmse = C_mmse.*signal_rx;% MMSE detection
            Dk_est_mmse(:,pilot_position) = [];% eliminate data on pilot position
            %---16QAM Demodulation---
            for l =1:(ncarriers-256)% 256 data eliminated
                C_Dk_zf = sqrt(abs(C-Dk_est_zf(l)).^2);% compute the Euclidean distance of every symbols to constellation points
                C_Dk_mmse = sqrt(abs(C-Dk_est_mmse(l)).^2);
                [a,index_zf] = min(C_Dk_zf);% map  points to constellation points with shortest distance
                [b,index_mmse] = min(C_Dk_mmse);
                rx_zf(l) = index_zf-1;% the original signal were 0-15
                rx_mmse(l) = index_mmse-1;
            end
            %===BER===
            tx(:,pilot_position) = [];% eliminate data on pilot position
            %---zf---
            bit_errors_zf = biterr(rx_zf,tx);% compute the Euclidean distance of every symbols to constellation points
            total_errors_zf = total_errors_zf+bit_errors_zf;% sum errors in average loops
            %---mmse---
            bit_errors_mmse = biterr(rx_mmse,tx);
            total_errors_mmse = total_errors_mmse+bit_errors_mmse;
            
            %===MSE===
            MSE_symbol(nav) = sum(abs(Hk-Hk_est).^2);
        end
        MSE(m,n) = 10.*log10((1./ncarriers).*sum(MSE_symbol)/Nav);
    end
    ber_sim_16qam_mp_zf(m) = total_errors_zf/(ncarriers*symbol_bits-256)/Nav;% compute BER
    ber_sim_16qam_mp_mmse(m) = total_errors_mmse/(ncarriers*symbol_bits-256)/Nav;
    
    toc
end

%===Plot results===
% Hk and estimated Hk
figure()
plot(abs(Hk),'r-');
hold on;
plot(abs(Hk_est),'b--');
title('SNR=30dB, K=64');
legend('Hk','Estimated Hk');
xlim([0 2048]);
grid on;
hold off;
% BER
figure()
semilogy(SNR_mp,Pb_MP,'b-',SNR_mp,ber_sim_16qam_mp_zf,'o-g',SNR_mp,ber_sim_16qam_mp_mmse,'*-y');
legend('16 QAM MP(Semi.)','16 QAM ZF MP(Sim.)','16 QAM MMSE MP(Sim.)');
xlabel('Eb/N0');
ylabel('BER')
xlim([0 30]);
ylim([1e-5 1]);
grid on;
% MSE
figure()
plot(SNR_mp,MSE(:,1),'b',SNR_mp,MSE(:,2),'g',SNR_mp,MSE(:,3),'r',SNR_mp,MSE(:,4),'c',SNR_mp,MSE(:,5),'m');
legend('K=0','K=8','K=16','K=32','K=64');
xlabel('Eb/N0(dB)');
ylabel('MSE(dB)');
xlim([0 30]);
ylim([-40 -10]);
grid on;
toc

%===Save results===
save MP_UnknownCH;