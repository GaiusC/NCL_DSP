clear, clc, close all
RN=sum(100*clock);
RS=RandStream('mt19937ar','seed',RN);
RandStream.setGlobalStream(RS); % obtain different random values

% Simulation parameter definitions
M=16;    % Constellation size
N=2048;  %number of subcarriers
C=[-3+3j,-3+1j,-3-3j,-3-1j,-1+3j,-1+1j,-1-3j,-1-1j,1+3j,1+1j,1-3j,1-1j,3+3j,3+1j,3-3j,3-1j];
SNR_dB=0:1:30;%Eb/N0
L=length(SNR_dB);
SNR_lin=10.^(SNR_dB/10);

Nt=2; % Number of Tx antennas
Nr=2; % Number of Rx antennas
H=(randn(Nr,Nt)+1j*randn(Nr,Nt))/sqrt(2);
Es=sum(abs(C).^2)/M;% average symbol energy
Eb_MI=Nt*Es/log2(M);
NO_MI=Eb_MI./SNR_lin;
sigma_MI=sqrt(NO_MI/2);
gamma_s = (Es.^2)./(sigma_MI.^2);
Nav=1000;
total_bits=Nav*log2(M);

% Equalizer zf
if Nt==Nr
    Czf=inv(H);
    
elseif  Nr>Nt
    Czf=pinv(H);
    %             Czf=inv(conj(H').*H).*conj(H');
end
% Channel matrix
% SNR loop

for m=1:L
    % Reset error count per SNR
    total_errors=0;
    total_errors2=0;
    Dest1=zeros(1,Nt);
    Dest2=zeros(1,Nt);

    % Averaging loop
    for n=1:Nav
        
        % Transmitter
        bk=randi([0 M-1],1,Nt);
        Dk=C(bk+1).';
        
        % Communications channel
        w=sigma_MI(m)*(randn(Nr,1)+1j*randn(Nr,1));
        
        % Receiver
        r=H*Dk+w;
        
        % Equalizer mmse
        Cmmse=inv((conj(H.'))*H+1/gamma_s(m))*(conj(H.'));
        
        Dest1=Czf*r;
        Dest2=Cmmse*r;
        % Bit error computation
        bk_est=zeros(1,Nt);
        bk_est2=zeros(1,Nt);
        
        for p=1:Nt
            [val,pos]=min(abs(C-Dest1(p)));
            bk_est(p)=pos-1;
            
            [val2,pos2]=min(abs(C-Dest2(p)));
            bk_est2(p)=pos2-1;
        end
        
        new_errors=biterr(bk,bk_est);
        new_errors2=biterr(bk,bk_est2);
        
        total_errors=total_errors+new_errors;
        total_errors2=total_errors2+new_errors2;
    end
    % Compute BER
    BER_zf(m)=total_errors/total_bits;
    BER_mmse(m)=total_errors2/total_bits;
    
end
%scatterplot(D);
% Store results in
semilogy(SNR_dB,BER_zf,'r-',SNR_dB,BER_mmse,'b-');
xlabel('E_b/N_0 (dB)')
ylabel('BER')
legend('zf','mmse')

axis([0 30 1e-5 1]);