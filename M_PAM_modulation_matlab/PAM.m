function [output_signal,input_symbols,output_symbols] = PAM(input_signal,M,Lb,Tc,Tsample,Tsymbol,SNR,coding)

% This function performs an M-PAM modulation in a sequence of sending bits
% through a channel.
%
% The function's inputs are :
% input_signal : the input sequence of bits to be transmitted
% M : the order of the M-PAM ( M = 4 or M = 8 )
% Lb : the number of the transmitting bits
% Tc : the time period of the carry
% Tsample : the time period of sampling
% Tsymbol : the time period of symbol
% SNR : the Signal-to-Noise-Ratio ( in db )
% coding : the type of bits' coding ( choose 'gray' or 'normal' )
%
% The function's outputs are :
% output_signal : the output sequence of bits
% input_symbols : the symbols' sequence before being transmitted
% output_symbols : the symbols' sequence after being transmitted

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
    
    % Normal coding for 4-PAM
    if(M==4 && coding == 'normal') 
        if(sum(sym)==0)
            Sm(i,1)=-3;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=-1;
            continue;
        elseif(sum(sym)==2)
            Sm(i,1)=3;
            continue;
        else
            Sm(i,1)=1;
            continue;
        end
    end
    
    % Gray coding for 4-PAM
    if(M==4 && coding == 'gray') 
        if(sum(sym)==0)
            Sm(i,1)=-3;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=-1;
            continue;
        elseif(sum(sym)==2)
            Sm(i,1)=1;
            continue;
        else
            Sm(i,1)=3;
            continue;
        end
    end
    
    % Normal Coding for 8-PAM
    if(M==8 && coding == 'normal') 
        if(sum(sym)==0)
            Sm(i,1)=-7;
            continue;
        elseif(sum(sym)==1 && sym(3,1)==1)
            Sm(i,1)=-5;
            continue;
        elseif(sum(sym)==2 && sym(1,1)==0)
            Sm(i,1)=-1;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=-3;
            continue;
        elseif(sum(sym)==2 && sym(3,1)==0)
            Sm(i,1)=5;
            continue;
        elseif(sum(sym)==3)
            Sm(i,1)=7;
            continue;
        elseif(sum(sym)==2 && sym(2,1)==0)
            Sm(i,1)=3;
            continue;
        else
            Sm(i,1)=1;
            continue;
        end
    end
    
    % Gray Coding for 8-PAM
    if(M==8 && coding == 'gray') 
        if(sum(sym)==0)
            Sm(i,1)=-7;
            continue;
        elseif(sum(sym)==1 && sym(3,1)==1)
            Sm(i,1)=-5;
            continue;
        elseif(sum(sym)==2 && sym(1,1)==0)
            Sm(i,1)=-3;
            continue;
        elseif(sum(sym)==1 && sym(2,1)==1)
            Sm(i,1)=-1;
            continue;
        elseif(sum(sym)==2 && sym(3,1)==0)
            Sm(i,1)=1;
            continue;
        elseif(sum(sym)==3)
            Sm(i,1)=3;
            continue;
        elseif(sum(sym)==2 && sym(2,1)==0)
            Sm(i,1)=5;
            continue;
        else
            Sm(i,1)=7;
            continue;
        end
    end
    
end
input_symbols = Sm;

% PAM modulator
sm=zeros(size*Tsymbol,1);
for i=1:size
    for j=1:Tsymbol
        sm(j+Tsymbol*(i-1),1)=sqrt(2/Tsymbol)*cos(2*pi*(1/Tc)*j)*Sm(i,1);
    end
end

% AWGN channel
Eb = 1/(log2(M));
No = Eb/(10^(SNR/10));
s2 = No/2;
noise = sqrt(s2)*randn(size*Tsymbol,1);
sm = sm+noise;

% PAM demodulator
Sm = zeros(size*Tsymbol,1);
for i=1:size
    for j=1:Tsymbol
        Sm(j+Tsymbol*(i-1),1)=sqrt(2/Tsymbol)*cos(2*pi*(1/Tc)*j)*sm(j+Tsymbol*(i-1),1);
    end
end

% Demodulator's summer
Sm = reshape(Sm,[Tsymbol size]);
sm = sum(Sm);
sm=sm';

% Matching of symbols
if M==4 
    symbols = [-3 -1 1 3];
else
    symbols = [-7 -5 -3 -1 1 3 5 7];
end
Sm = zeros(size,1);
distances = zeros(size,length(symbols));
for i =1 : size
    for j =1 : M
        distances(i,j)=abs(sm(i,1)-symbols(1,j));
    end
    [~,index_min]=min(distances(i,1:M));
    Sm(i,1) = symbols(1,index_min); 
end
output_symbols = Sm;

% Demapper
if M==4 && coding=='normal'
    output_signal=[];
    for i=1:size
        if Sm(i,1)==-3
            output_signal = [output_signal;0;0];
        elseif Sm(i,1)==-1
            output_signal = [output_signal;0;1];
        elseif Sm(i,1)==1
            output_signal = [output_signal;1;0];
        else
            output_signal = [output_signal;1;1];
        end
    end
end

if M==4 && coding=='gray'
    output_signal=[];
    for i=1:size
        if Sm(i,1)==-3
            output_signal = [output_signal;0;0];
        elseif Sm(i,1)==-1
            output_signal = [output_signal;0;1];
        elseif Sm(i,1)==1
            output_signal = [output_signal;1;1];
        else
            output_signal = [output_signal;1;0];
        end
    end
end

if M==8 && coding=='normal'
    output_signal=[];
    for i=1:size
        if Sm(i,1)==-7
            output_signal = [output_signal;0;0;0];
        elseif Sm(i,1)==-5
            output_signal = [output_signal;0;0;1];
        elseif Sm(i,1)==-3
            output_signal = [output_signal;0;1;0];
        elseif Sm(i,1)==-1
            output_signal = [output_signal;0;1;1];
        elseif Sm(i,1)==1
            output_signal = [output_signal;1;0;0];
        elseif Sm(i,1)==3
            output_signal = [output_signal;1;0;1];
        elseif Sm(i,1)==5
            output_signal = [output_signal;1;1;0];
        else
            output_signal = [output_signal;1;1;1];
        end
    end
end

if M==8 && coding=='gray'
    output_signal=[];
    for i=1:size
        if Sm(i,1)==-7
            output_signal = [output_signal;0;0;0];
        elseif Sm(i,1)==-5
            output_signal = [output_signal;0;0;1];
        elseif Sm(i,1)==-3
            output_signal = [output_signal;0;1;1];
        elseif Sm(i,1)==-1
            output_signal = [output_signal;0;1;0];
        elseif Sm(i,1)==1
            output_signal = [output_signal;1;1;0];
        elseif Sm(i,1)==3
            output_signal = [output_signal;1;1;1];
        elseif Sm(i,1)==5
            output_signal = [output_signal;1;0;1];
        else
            output_signal = [output_signal;1;0;0];
        end
    end
end

end
