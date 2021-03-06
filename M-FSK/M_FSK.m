function [output_signal,input_signal,W] = M_FSK(input_signal,M,Lb,Tc,Tsample,Tsymbol,SNR,Es,coding)

% This function performs an M-FSK modulation in a sequence of sending bits
% through a channel.
%
% The function's inputs are :
% input_signal : the input sequence of bits to be transmitted
% M : the order of the M_FSK ( M = 4 or M = 8 )
% Lb : the number of the transmitting bits
% Tc : the time period of the carry
% Tsample : the time period of sampling
% Tsymbol : the time period of symbol
% SNR : the Signal-to-Noise-Ratio ( in db )
% Es : the energy of each symbol ( we set it equal to 1, normalization )
% coding : the type of bits' coding ( choose 0 for normal code or 1 for gray code )
%
% The function's outputs are :
% output_signal : the output sequence of bits
% input_signal : the input sequence of bits with the extra 0s at the end

% Adding 0s to the end of input's bit sequence
count = 0;
while(true)
   if(mod(length(input_signal(:,1)),log2(M))==0)
       break;
   end
   input_signal=[input_signal;0];
   count=count+1;
end
Lb = Lb+count;

% Periods' normalization
Tsample_norm=1;
k = Tsample_norm/Tsample;
Tc = round(k*Tc);
Tsymbol = round(k*Tsymbol);

% Number of symbols' sequence
size = Lb/(log2(M));

% Mapper ( matching the bits to symbols )
Sm = zeros(size,1);
c1=1;
c2=0;
for i=1:(size)
    sym = zeros(log2(M),1);
    for j=c1:(c1+log2(M)-1)
        sym (j-c2*log2(M),1)=input_signal(j,1);
        c1=c1+1;
    end
    c2=c2+1;
    
    % coding for BFSK
    if(M==2 && coding == 0)
       if(sym == 0)
           Sm(i,1)=-1;
       else
           Sm(i,1)=1;
       end
    end
    
    % Normal coding for 4-PSK
    if(M==4 && coding == 0) 
        if(sum(sym)==0)
            Sm(i,1)=1;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=2;
            continue;
        elseif(sum(sym)==2)
            Sm(i,1)=4;
            continue;
        else
            Sm(i,1)=3;
            continue;
        end
    end
    
    % Gray coding for 4-PSK
    if(M==4 && coding == 1) 
        if(sum(sym)==0)
            Sm(i,1)=0;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=1;
            continue;
        elseif(sum(sym)==2)
            Sm(i,1)=2;
            continue;
        else
            Sm(i,1)=3;
            continue;
        end
    end
    
    % Normal Coding for 8-PSK
    if(M==8 && coding == 0) 
        if(sum(sym)==0)
            Sm(i,1)=0;
            continue;
        elseif(sum(sym)==1 && sym(3,1)==1)
            Sm(i,1)=1;
            continue;
        elseif(sum(sym)==2 && sym(1,1)==0)
            Sm(i,1)=3;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=2;
            continue;
        elseif(sum(sym)==2 && sym(3,1)==0)
            Sm(i,1)=6;
            continue;
        elseif(sum(sym)==3)
            Sm(i,1)=7;
            continue;
        elseif(sum(sym)==2 && sym(2,1)==0)
            Sm(i,1)=5;
            continue;
        else
            Sm(i,1)=4;
            continue;
        end
    end
    
    % Gray Coding for 8-PSK
    if(M==8 && coding == 1) 
        if(sum(sym)==0)
            Sm(i,1)=0;
            continue;
        elseif(sum(sym)==1 && sym(3,1)==1)
            Sm(i,1)=1;
            continue;
        elseif(sum(sym)==2 && sym(1,1)==0)
            Sm(i,1)=2;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=3;
            continue;
        elseif(sum(sym)==2 && sym(3,1)==0)
            Sm(i,1)=4;
            continue;
        elseif(sum(sym)==3)
            Sm(i,1)=5;
            continue;
        elseif(sum(sym)==2 && sym(2,1)==0)
            Sm(i,1)=6;
            continue;
        else
            Sm(i,1)=7;
            continue;
        end
    end
    
    % Normal Coding for 16-PSK
    if(M==16 && coding==0)
        if(sum(sym)==0)
            Sm(i,1)=0;
        elseif(sum(sym)==1 && sym(4,1)==1)
            Sm(i,1)=1;
        elseif(sum(sym)==1 && sym(3,1)==1)
            Sm(i,1)=2;
        elseif(sum(sym)==2 && sym(4,1)==1 && sym(3,1)==1)
            Sm(i,1)=3; 
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=4;
        elseif(sum(sym)==2 && sym(4,1)==1 && sym(2,1)==1)
            Sm(i,1)=5;
        elseif(sum(sym)==2 && sym(2,1)==1 && sym(3,1)==1)
            Sm(i,1)=6;
        elseif(sum(sym)==3 && sym(1,1)==0)
            Sm(i,1)=7;
        elseif(sum(sym)==1 && sym(1,1)==1)
            Sm(i,1)=8;
        elseif(sum(sym)==2 && sym(1,1)==1 && sym(4,1)==1)
            Sm(i,1)=9;
        elseif(sum(sym)==2 && sym(1,1)==1 && sym(3,1)==1)
            Sm(i,1)=10;
        elseif(sum(sym)==3 && sym(2,1)==0)
            Sm(i,1)=11;
        elseif(sum(sym)==2 && sym(3,1)==1 && sym(4,1)==1)
            Sm(i,1)=12;
        elseif(sum(sym)==3 && sym(3,1)==0)
            Sm(i,1)=13;
        elseif(sum(sym)==3 && sym(4,1)==0)
            Sm(i,1)=14;
        else
            Sm(i,1)=15;
        end
    end
    
    % Gray Coding for 16-PSK
    if(M==16 && coding==1)
        if(sum(sym)==0)
            Sm(i,1)=0;
        elseif(sum(sym)==1 && sym(4,1)==1)
            Sm(i,1)=1;
        elseif(sum(sym)==1 && sym(3,1)==1)
            Sm(i,1)=3;
        elseif(sum(sym)==2 && sym(4,1)==1 && sym(3,1)==1)
            Sm(i,1)=2; 
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=6;
        elseif(sum(sym)==2 && sym(4,1)==1 && sym(2,1)==1)
            Sm(i,1)=7;
        elseif(sum(sym)==2 && sym(2,1)==1 && sym(3,1)==1)
            Sm(i,1)=5;
        elseif(sum(sym)==3 && sym(1,1)==0)
            Sm(i,1)=4;
        elseif(sum(sym)==1 && sym(1,1)==1)
            Sm(i,1)=12;
        elseif(sum(sym)==2 && sym(1,1)==1 && sym(4,1)==1)
            Sm(i,1)=13;
        elseif(sum(sym)==2 && sym(1,1)==1 && sym(3,1)==1)
            Sm(i,1)=15;
        elseif(sum(sym)==3 && sym(2,1)==0)
            Sm(i,1)=14;
        elseif(sum(sym)==2 && sym(3,1)==1 && sym(4,1)==1)
            Sm(i,1)=10;
        elseif(sum(sym)==3 && sym(3,1)==0)
            Sm(i,1)=11;
        elseif(sum(sym)==3 && sym(4,1)==0)
            Sm(i,1)=9;
        else
            Sm(i,1)=8;
        end
    end    
end 


% % % FSK modulator % % %
% Carrier frequency modulation
sm=zeros(size*Tsymbol,M);
for i = 1:size
    for j = 1:Tsymbol
        for z = 1:M
            sm(j+Tsymbol*(i-1),z)=Sm(i,1)*sqrt((2*Es)/Tsymbol)*cos(2*pi*((1/Tc)+((z-1)/Tsymbol))*j);
        end
    end
end

Sm = zeros(size*Tsymbol,1);
for i = 1:M
    Sm = Sm + sm(:,i);
end

% AWGN channel
Eb = 1/(log2(M));
No = Eb/(10^(SNR/10));
s2 = No/2; % Variance of additive white gaussian noise
noise = sqrt(s2)*randn(size*Tsymbol,1);
Sm = Sm+noise;

% % % FSK demodulator % % %
sm = zeros(size*Tsymbol,M);
for i = 1:size
    for j = 1:Tsymbol
        for z = 1:M
            sm(j+Tsymbol*(i-1),z)=sqrt((2*Es)/Tsymbol)*Sm(j+Tsymbol*(i-1),1)*cos(2*pi*((1/Tc)+((z-1)/Tsymbol))*j);
        end
    end
end

% Demodulator's summer
Sm = zeros(size,M);
for z = 1:M
    for i = 1:size
        for j = 1:Tsymbol
            Sm(i,z)=Sm(i,z)+sm(j+Tsymbol*(i-1),z);
        end
    end
end
W = Sm;
% % % Matching of symbols % % %

% Matching
sm = zeros(size,1);
for i =1 : size
    [~,index_max]=max(Sm(i,:));
    sm(i,1) = index_max;
end
output_symbols = sm;
Sm = sm;

% % % Demapper % % %
% Coding for M=2
if M==2
    output_signal=[];
    for i=1:size
        if Sm(i,1)==1
            output_signal = [output_signal;0];
        else
            output_signal = [output_signal;1];
        end
    end
end
% Normal coding for M=4
if M==4 && coding==0
    output_signal=[];
    for i=1:size
        if Sm(i,1)==1
            output_signal = [output_signal;0;0];
        elseif Sm(i,1)==2
            output_signal = [output_signal;0;1];
        elseif Sm(i,1)==3
            output_signal = [output_signal;1;0];
        else
            output_signal = [output_signal;1;1];
        end
    end
end
% Gray coding for M=4
if M==4 && coding==1
    output_signal=[];
    for i=1:size
        if Sm(i,1)==0
            output_signal = [output_signal;0;0];
        elseif Sm(i,1)==1
            output_signal = [output_signal;0;1];
        elseif Sm(i,1)==2
            output_signal = [output_signal;1;1];
        else
            output_signal = [output_signal;1;0];
        end
    end
end
% Normal coding for M=8
if M==8 && coding==0
    output_signal=[];
    for i=1:size
        if Sm(i,1)==0
            output_signal = [output_signal;0;0;0];
        elseif Sm(i,1)==1
            output_signal = [output_signal;0;0;1];
        elseif Sm(i,1)==2
            output_signal = [output_signal;0;1;0];
        elseif Sm(i,1)==3
            output_signal = [output_signal;0;1;1];
        elseif Sm(i,1)==4
            output_signal = [output_signal;1;0;0];
        elseif Sm(i,1)==5
            output_signal = [output_signal;1;0;1];
        elseif Sm(i,1)==6
            output_signal = [output_signal;1;1;0];
        else
            output_signal = [output_signal;1;1;1];
        end
    end
end
% Gray coding for M=8
if M==8 && coding==1
    output_signal=[];
    for i=1:size
        if Sm(i,1)==0
            output_signal = [output_signal;0;0;0];
        elseif Sm(i,1)==1
            output_signal = [output_signal;0;0;1];
        elseif Sm(i,1)==2
            output_signal = [output_signal;0;1;1];
        elseif Sm(i,1)==3
            output_signal = [output_signal;0;1;0];
        elseif Sm(i,1)==4
            output_signal = [output_signal;1;1;0];
        elseif Sm(i,1)==5
            output_signal = [output_signal;1;1;1];
        elseif Sm(i,1)==6
            output_signal = [output_signal;1;0;1];
        else
            output_signal = [output_signal;1;0;0];
        end
    end
end
% Normal coding for M=16
if M == 16 && coding == 0
    output_signal=[];
    for i=1:size
        if Sm(i,1)==0
            output_signal = [output_signal;0;0;0;0];
        elseif Sm(i,1)==1
            output_signal = [output_signal;0;0;0;1];
        elseif Sm(i,1)==2
            output_signal = [output_signal;0;0;1;0];
        elseif Sm(i,1)==3
            output_signal = [output_signal;0;0;1;1];
        elseif Sm(i,1)==4
            output_signal = [output_signal;0;1;0;0];
        elseif Sm(i,1)==5
            output_signal = [output_signal;0;1;0;1];
        elseif Sm(i,1)==6
            output_signal = [output_signal;0;1;1;0];
        elseif Sm(i,1)==7
            output_signal = [output_signal;0;1;1;1];
        elseif Sm(i,1)==8
            output_signal = [output_signal;1;0;0;0];
        elseif Sm(i,1)==9
            output_signal = [output_signal;1;0;0;1];
        elseif Sm(i,1)==10
            output_signal = [output_signal;1;0;1;0];
        elseif Sm(i,1)==11
            output_signal = [output_signal;1;0;1;1];
        elseif Sm(i,1)==12
            output_signal = [output_signal;1;1;0;0];
        elseif Sm(i,1)==13
            output_signal = [output_signal;1;1;0;1];
        elseif Sm(i,1)==14
            output_signal = [output_signal;1;1;1;0];
        else
            output_signal = [output_signal;1;1;1;1];
        end
    end
end
% Gray coding for M=16
if M == 16 && coding == 1
    output_signal=[];
    for i=1:size
        if Sm(i,1)==0
            output_signal = [output_signal;0;0;0;0];
        elseif Sm(i,1)==1
            output_signal = [output_signal;0;0;0;1];
        elseif Sm(i,1)==3
            output_signal = [output_signal;0;0;1;0];
        elseif Sm(i,1)==2
            output_signal = [output_signal;0;0;1;1];
        elseif Sm(i,1)==6
            output_signal = [output_signal;0;1;0;0];
        elseif Sm(i,1)==7
            output_signal = [output_signal;0;1;0;1];
        elseif Sm(i,1)==5
            output_signal = [output_signal;0;1;1;0];
        elseif Sm(i,1)==4
            output_signal = [output_signal;0;1;1;1];
        elseif Sm(i,1)==12
            output_signal = [output_signal;1;0;0;0];
        elseif Sm(i,1)==13
            output_signal = [output_signal;1;0;0;1];
        elseif Sm(i,1)==15
            output_signal = [output_signal;1;0;1;0];
        elseif Sm(i,1)==14
            output_signal = [output_signal;1;0;1;1];
        elseif Sm(i,1)==10
            output_signal = [output_signal;1;1;0;0];
        elseif Sm(i,1)==11
            output_signal = [output_signal;1;1;0;1];
        elseif Sm(i,1)==9
            output_signal = [output_signal;1;1;1;0];
        else
            output_signal = [output_signal;1;1;1;1];
        end
    end
end
end