clc;clear all;close all;
nsymbol=2048;%表示一共有多少个符号，这里定义100000个符号
M=16;%M表示QAM调制的阶数,表示16QAM，16QAM采用格雷映射(所有星座点图均采用格雷映射)
N=64;
graycode=[0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];%格雷映射编码规则
graycode1=[0 1 3 2 6 7 5 4 8 9 11 10 14 15 13 12 24 25 27 26 30 31 29 28 16 17 19 18 22 23 21 20 48 49 51 50 54 55 53 52 56 57 59 58 62 63 61 60 40 41 43 42 46 47 45 44 32 33 35 34 38 39 37 36];%格雷映射十进制的表示
EsN0=5:20;%信噪比范围
snr1=10.^(EsN0/10);%将db转换为线性值
msg=randi([0,M-1],1,nsymbol);%0到15之间随机产生一个数,数的个数为：1乘nsymbol，得到原始数据
msg1=graycode(msg+1);%对数据进行格雷映射
msgmod=qammod(msg1,M);%调用matlab中的qammod函数，16QAM调制方式的调用(输入0到15的数，M表示QAM调制的阶数)得到调制后符号

scatterplot(msgmod);%调用matlab中的scatterplot函数,画星座点图
spow=norm(msgmod).^2/nsymbol;%取a+bj的模.^2得到功率除整个符号得到每个符号的平均功率
%64QAM
nsg=randi([0,N-1],1,nsymbol);
nsg1=graycode1(nsg+1);
nsgmod=qammod(nsg1,N);
scatterplot(nsgmod);%调用matlab中的scatterplot函数,画星座点图
spow1=norm(nsgmod).^2/nsymbol;

for i=1:length(EsN0)
    sigma=sqrt(spow/(2*snr1(i)));%16QAM根据符号功率求出噪声的功率
    sigma1=sqrt(spow1/(2*snr1(i)));%64QAM根据符号功率求出噪声的功率
    rx=msgmod+sigma*(randn(1,length(msgmod))+1i*randn(1,length(msgmod)));%16QAM混入高斯加性白噪声
    rx1=nsgmod+sigma*(randn(1,length(nsgmod))+1i*randn(1,length(nsgmod)));%64QAM混入高斯加性白噪声
    y=qamdemod(rx,M);%16QAM的解调
   y1=qamdemod(rx1,N);%64QAM的解调
   decmsg=graycode(y+1);%16QAM接收端格雷逆映射，返回译码出来的信息，十进制
   decnsg=graycode1(y1+1);%64QAM接收端格雷逆映射
   [err1,ber(i)]=biterr(msg,decmsg,log2(M));%g一个符号四个比特，比较发送端信号msg和解调信号decmsg转换为二进制，ber(i)错误的比特率
   [err2,ser(i)]=symerr(msg,decmsg);%16QAM求实际误码率
   [err1,ber1(i)]=biterr(nsg,decnsg,log2(N));
   [err2,ser1(i)]=symerr(nsg,decnsg);%64QAM求实际误码率
end
%16QAM
scatterplot(rx);%调用matlab中的scatterplot函数,画rx星座点图
p = 2*(1-1/sqrt(M))*qfunc(sqrt(3*snr1/(M-1)));
ser_theory=1-(1-p).^2;%16QAM理论误码率
ber_theory=1/log2(M)*ser_theory;

%64QAM
scatterplot(rx1);
p1=2*(1-1/sqrt(N))*qfunc(sqrt(3*snr1/(N-1)));
ser1_theory=1-(1-p1).^2;%64QAM理论误码率
ber1_theory=1/log2(N)*ser1_theory;%得到误比特率

%绘图
figure()
semilogy(EsN0,ber,"o", EsN0, ser, "*",EsN0, ser_theory, "-", EsN0, ber_theory, "-");
title("16-QAM载波调制信号在AWGN信道下的误比特率性能")
xlabel("EsN0");
ylabel("误比特率和误符号率");
legend("误比特率", "误符号率","理论误符号率","理论误比特率");
%阶数不同,16和64QAM调制信号在AWGN信道的性能比较
figure()
semilogy(EsN0,ser_theory,'o',EsN0,ser1_theory,'o');%ber ser比特仿真值 ser1理论误码率 ber1理论误比特率
title('16和64QAM调制信号在AWGN信道的性能比较');grid;
xlabel('Es/N0(dB)');%性躁比
ylabel('误码率');%误码率
legend('16QAM理论误码率','64QAM理论误码率');