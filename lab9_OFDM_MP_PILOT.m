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
Eb=Es/symbol_bits;
N0=Eb./SNR_mp_lin;
sigma=sqrt(N0/2);
gama_s=(Es.^2)./(sigma.^2);
pilot_content = 3+3i;
pilot_position = 1:8:2048;

%===Channel Energy Normalization===
P=[0 -0.9 -4.9 -8 -7.8 -23.9];
h_l=sqrt(10.^(P/10));
U=sqrt(sum(abs(h_l).^2));
hl=h_l./U;
H=[hl(1),zeros(1,1),hl(2),zeros(1,6),hl(3),zeros(1,4),hl(4),zeros(1,11),hl(5),zeros(1,15),hl(6)];
Hk=fft(H,ncarriers);


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
Nav = 50;
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
            tx = randi([0 M-1],1,ncarriers);
            Dn_16qam_mod = C(tx+1);
            %---Insert pilot---
            Dn_pilot = Dn_16qam_mod;
            Dn_pilot(1,pilot_position) = pilot_content;
            %---ifft---
            signal_tx = ifft(Dn_pilot)*sqrt(ncarriers);
            %---cp---
            xn = [signal_tx(1,ncarriers-cp_length(n)+1:end),signal_tx];
            
            %===Channel===
            [yn,zf] = filter(H,1,xn,zf);
            awgn_noise = sigma(m)*(randn(1,ncarriers+cp_length(n))+1i*randn(1,ncarriers+cp_length(n)));
            signal_awgn = yn+awgn_noise;
            
            %===Rx===
            %---Remove Cyclic Prefix---
            rx_cp = signal_awgn(:,cp_length(n)+1:end);
            %---fft---
            signal_rx = fft(rx_cp)/sqrt(ncarriers);
            %---equalizer---
            Hp_est = signal_rx(pilot_position)/pilot_content;
            k = 0;
            for pp = 1:255
                for ii = 0:7
                    k = k+1;
                    Hk_est(k) = Hp_est(pp)+(ii/8)*(Hp_est(pp+1)-Hp_est(pp));
                end
            end
            for iii = 8:15
                k = k+1;
                Hk_est(k) = Hp_est(pp)+(iii/8)*(Hp_est(pp+1)-Hp_est(pp));
            end
            C_zf = conj(Hk_est)./(abs(Hk_est).^2);
            C_mmse = conj(Hk_est)./((abs(Hk_est).^2)+(1./gama_s(m)));
            Dk_est_zf = C_zf.*signal_rx;
            Dk_est_zf(:,pilot_position) = [];
            Dk_est_mmse = C_mmse.*signal_rx;
            Dk_est_mmse(:,pilot_position) = [];
            %---16QAM Demodulation---
            for l =1:(ncarriers-256)
                C_Dk_zf = sqrt(abs(C-Dk_est_zf(l)).^2);
                C_Dk_mmse = sqrt(abs(C-Dk_est_mmse(l)).^2);
                [a,index_zf] = min(C_Dk_zf);
                [b,index_mmse] = min(C_Dk_mmse);
                rx_zf(l) = index_zf-1;
                rx_mmse(l) = index_mmse-1;
            end
            %===BER===
            tx(:,pilot_position) = [];
            %---zf---
            bit_errors_zf = biterr(rx_zf,tx);
            total_errors_zf = total_errors_zf+bit_errors_zf;
            %---mmse---
            bit_errors_mmse = biterr(rx_mmse,tx);
            total_errors_mmse = total_errors_mmse+bit_errors_mmse;
            
            %===MSE===
            MSE_symbol(nav) = sum(abs(Hk-Hk_est).^2);
        end
        MSE(m,n) = 10.*log10((1./ncarriers).*sum(MSE_symbol)/Nav);
    end
    ber_sim_16qam_mp_zf(m) = total_errors_zf/(ncarriers*symbol_bits-256)/Nav;
    ber_sim_16qam_mp_mmse(m) = total_errors_mmse/(ncarriers*symbol_bits-256)/Nav;
    
    toc
end

%===Plot results===
% Czf and Cmmes
figure()
plot(Dk_est_zf,'o');
title('Czf');
figure()
plot(Dk_est_mmse,'o');
title('Cmmse');
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

%===Save results===
save lab9_OFDM_MP_PILOT;
toc