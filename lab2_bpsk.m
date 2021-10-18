clear, clc, close all
tic
%--BPSK in AWGN--
gamma_dB=-5:1:10;
gamma_lin=10.^(gamma_dB/10);
ber_AWGN=qfunc(sqrt(2*gamma_lin));

%--BPSK in Rayleigh--
Gamma_dB=-5:5:50;
Gamma_lin=10.^(Gamma_dB/10);
ber_Rayleigh=(1-sqrt(Gamma_lin./(Gamma_lin+1)))/2;

%--Semi-analytical--
N=1e6;
gamma_dB_2=-5:5:50;
gamma_lin=10.^(gamma_dB_2/10);
L=length(gamma_dB_2);
ber_Semi_Rayleigh=zeros(1,L);
for i=1:L
    h=1/sqrt(2)*(randn(1,N)+1j*randn(1,N));
    Pe=qfunc(sqrt(2*gamma_lin(i)*abs(h).^2));
    ber_Semi_Rayleigh(i)=sum(Pe)/N;
    toc
end

%Full Simulation
nsymbol=1e6;
gamma_dB_3=-5:5:45;
gamma_lin=10.^(gamma_dB_3/10);
C=[1 -1];
ber_full=zeros(length(gamma_dB_3),1);
dn_r=zeros(1,nsymbol);
bn_r=zeros(1,nsymbol);

for m=1:length(gamma_dB_3)
    
    total_errors=0;
    
    %Tx
    bn=randi([0 1],1,nsymbol);
    dn=C(bn+1);
    
    %Channel
    hn=1/sqrt(2)*(randn(1,nsymbol)+1i*randn(1,nsymbol));
    sigma=sqrt(0.5/gamma_lin(m));
    wn=sigma*(randn(1,nsymbol)+1i*randn(1,nsymbol));
    rn=hn.*dn+wn;
    %Rx
    for n=1:nsymbol
        dn_r(n)=real(rn(n)./hn(n));
        if dn_r(n)<=0
            bn_r(n)=1;
        else
            bn_r(n)=0;
        end
    end
    %BER
    new_errors=biterr(bn_r,bn);
    total_errors=total_errors+new_errors;
    ber_full(m)=total_errors/nsymbol;
    toc
end

%Plot results(full)
semilogy(gamma_dB,ber_AWGN,'r-',Gamma_dB,ber_Rayleigh,'b--',gamma_dB_2,ber_Semi_Rayleigh,'o--',gamma_dB_3,ber_full,'go--');
xlabel("Average SNR ");
ylabel("Pe");
ylim([1e-6 1]);
grid on;
legend("BPSK in AWGN","BPSK in Rayleigh","Semi-analytial","Full Simulation");
toc
%save results
%save lab2_bpsk_results;
