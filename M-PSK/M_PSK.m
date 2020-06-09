function [output_signal,input_signal] = M_PSK(input_signal,M,Lb,Tc,Tsample,Tsymbol,SNR,Es,coding)

% This function performs an M-PSK modulation in a sequence of sending bits
% through a channel.
%
% The function's inputs are :
% input_signal : the input sequence of bits to be transmitted
% M : the order of the M_PSK ( M = 4 or M = 8 )
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
           Sm(i,1)=0;
       else
           Sm(i,1)=1;
       end
    end
    
    % Normal coding for 4-PSK
    if(M==4 && coding == 0) 
        if(sum(sym)==0)
            Sm(i,1)=0;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=1;
            continue;
        elseif(sum(sym)==2)
            Sm(i,1)=3;
            continue;
        else
            Sm(i,1)=2;
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
input_symbols = Sm;

% % % PSK modulator % % %
% Symbols
Sm1 = zeros(size,1);
Sm2 = zeros(size,1);
for i =1:size
    Sm1(i,1)=sqrt(Es)*cos((2*pi*Sm(i,1))/4);
    Sm2(i,1)=sqrt(Es)*sin((2*pi*Sm(i,1))/4);
end
% Carrier frequency modulation
sm1=zeros(size*Tsymbol,1);
sm2=zeros(size*Tsymbol,1);
for i = 1:size
    for j = 1:Tsymbol
        sm1(j+Tsymbol*(i-1),1)=sqrt((2*Es)/Tsymbol)*Sm1(i,1)*cos(2*pi*(1/Tc)*j);
        sm2(j+Tsymbol*(i-1),1)=sqrt((2*Es)/Tsymbol)*Sm2(i,1)*sin(2*pi*(1/Tc)*j);
    end
end
Sm = sm1 + sm2;

% AWGN channel
Eb = 1/(log2(M));
No = Eb/(10^(SNR/10));
s2 = No/2; % Variance of additive white gaussian noise
noise = sqrt(s2)*randn(size*Tsymbol,1);
Sm = Sm+noise;

% % % PSK demodulator % % %
Sm1 = zeros(size*Tsymbol,1);
Sm2 = zeros(size*Tsymbol,1);
for i = 1:size
    for j = 1:Tsymbol
        Sm1(j+Tsymbol*(i-1),1)=sqrt((2*Es)/Tsymbol)*Sm(j+Tsymbol*(i-1),1)*cos(2*pi*(1/Tc)*j);
        Sm2(j+Tsymbol*(i-1),1)=sqrt((2*Es)/Tsymbol)*Sm(j+Tsymbol*(i-1),1)*sin(2*pi*(1/Tc)*j);
    end
end

% Demodulator's summer
Sm1 = reshape(Sm1,[Tsymbol size]);
Sm2 = reshape(Sm2,[Tsymbol size]);
sm1 = sum(Sm1);
sm2 = sum(Sm2);

% Final Vector
r = [sm1;sm2];

% % % Matching of symbols % % %
% Making of symbols vector
symbols1=zeros(1,M);
symbols2=zeros(1,M);
symbols1(1,1)=cos(0);
symbols2(1,1)=sin(0);
for i= 1:(M-1)
    symbols1(1,i+1)=sqrt(Es)*cos((2*pi*i)/4);
    symbols2(1,i+1)=sqrt(Es)*sin((2*pi*i)/4);
end
symbols=[symbols1;symbols2];
% Matching
Sm = zeros(size,1);
distances = zeros(size,length(symbols(1,:)));
for i =1 : size
    for j =1 : M
        distances(i,j)=norm(r(:,i)-symbols(:,j));
    end
    [~,index_min]=min(distances(i,:));
    Sm(i,1) = index_min-1;
end
output_symbols = Sm;

% % % Demapper % % %
% Coding for M=2
if M==2
    output_signal=[];
    for i=1:size
        if Sm(i,1)==0
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
        if Sm(i,1)==0
            output_signal = [output_signal;0;0];
        elseif Sm(i,1)==1
            output_signal = [output_signal;0;1];
        elseif Sm(i,1)==2
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