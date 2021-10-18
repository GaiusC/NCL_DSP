clear, clc, close all
% Simulation parameter definitions
M=16; % Constellation size
nsymbol=2048;%define 2048 symbols
C=[-3+3j -3+1j -3-1j -3-3j -1+3j -1+1j -1-1j -1-3j 1+3j 1+1j 1-1j 1-3j 3+3j 3+1j 3-1j 3-3j];% Graycode rule
C_order=[0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];
spow=mean(abs(C).^2);
SNR=5:20;%range of SNR
SNR_lin=10.^(SNR/10);%convert dB to linear value

% SNR loop
for m=1:length(SNR)
    % Reset error count per SNR
    total_errors=0;
    %sigma(m)=sqrt(5/SNR_lin(m));
    Bk=randi([0,M-1],1,nsymbol);
    Bk_bin=double(de2bi(Bk))-48;
    Dk=C(Bk+1);
    
    
    % Averaging loop
    for n=1:nsymbol
        % Transmitter
        Bk(m,n)=randi([0,M-1]);
        Dk(m,n)=C(Bk(m,n)+1);
        if real(Dk(m,n))>0
            b0=1;
        else
            b0=0;
        end
        if real(Dk(m,n))>-2
            if real(Dk(m,n))<2
                b1=1;
            else
                b1=0;
            end
        else
            b1=0;
        end
        if imag(Dk(m,n))<0
            b2=1;
        else
            b2=0;
        end
        if imag(Dk(m,n))>-2
            if imag(Dk(m,n))<2
                b3=1;
            else
                b3=0;
            end
        else
            b3=0;
        end
        Dk_bin(n,:)=[b0 b1 b2 b3];
        Dk_bin_str(n,:)=num2str((Dk_bin(n,:))');
        Dk_dec(m,n)=bin2dec(Dk_bin_str(n,:));
        
        
        % Communications channel
        pow_Dk_bin(m,n)=norm(Dk(m,n)).^2;
        sigma(m,n)=sqrt(pow_Dk_bin(m,n)/(2*SNR_lin(m)));
        awgn_noise(m,n)=sigma(m,n)*(randn+1i*randn);
        
        % Receiver
        msg(m,n)=Dk(m,n)+awgn_noise(m,n);
        %msg(m,n)=awgn(Dk(m,n),SNR(m),'measured');
        if real(msg(m,n))>0
            b0=1;
        else
            b0=0;
        end
        if real(msg(m,n))>-2
            if real(msg(m,n))<2
                b1=1;
            else
                b1=0;
            end
        else
            b1=0;
        end
        if imag(msg(m,n))<0
            b2=1;
        else
            b2=0;
        end
        if imag(msg(m,n))>-2
            if imag(msg(m,n))<2
                b3=1;
            else
                b3=0;
            end
        else
            b3=0;
        end
        msg_bin(n,:)=[b0 b1 b2 b3];
        msg_bin_str(n,:)=num2str(msg_bin(n,:));
        msg_dec_col(n,:)=bin2dec(msg_bin_str(n,:));
        
        % Bit error computation
        new_errors=sum(Dk_bin(n,:)~=msg_bin(n,:));
        total_errors=total_errors+new_errors;
    end
    msg_dec(m,:)=msg_dec_col';
    demsg(m,:)=C_order(msg_dec(m,:)+1);
    
    % Compute BER
    BER(m)=total_errors/(log2(M)*nsymbol);
    BER1(m)=biterr(Bk(m,:),demsg(m,:),log2(M))/nsymbol/4;
end
%plot results
msg_025=Dk(randi([1 16]),:)+0.25*(randn(1,length(Dk))+1j*randn(1,length(Dk)));
scatterplot(Dk(1,:));
scatterplot(msg_025);

p = 2*(1-1/sqrt(M))*qfunc(sqrt(3*SNR_lin/(M-1)));
ser_theory=1-(1-p).^2;
ber_theory=1/log2(M)*ser_theory;
figure()
semilogy(SNR,BER1,"o", SNR, ber_theory, "-");
title("BER of 16-QAM in AWGN")
xlabel("EsN0");
ylabel("BER");
legend("BER", "Theoretical BER");

%save results
%save lab3_16qam_results
