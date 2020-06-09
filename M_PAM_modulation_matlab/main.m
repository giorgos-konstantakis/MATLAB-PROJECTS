%% 
clc; clear all;
%% Setting the values for the M_PAM 

% Number of bits
Lb = 100002;

% Input sequence
Binary_Input = randsrc(Lb,1,[0 1]);

% M_PAM's order
M=4; % choose  M=4(one symbol per 2 bits) or M=8(one symbol per 3 bits)

% The time periods are normalized , with setting Tsample=1
% Time period of symbol
Tsymbol=40;
% Time period of cosine's carrier
Tc=4;
% Time period of sampling
Tsample=1;

% Coding bits
coding = 'gray'

% SNR in db ( we set many SNR values , to calculate the BER diagram )
SNR = 0:2:20; %db

%% Extract the results of the M_PAM algorithm

% M-PAM application
[Output_Sequence_of_Bits_1,Input_Sequence_of_Symbols_1,Output_Sequence_of_Symbols_1] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(1));
[Output_Sequence_of_Bits_2,Input_Sequence_of_Symbols_2,Output_Sequence_of_Symbols_2] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(2));
[Output_Sequence_of_Bits_3,Input_Sequence_of_Symbols_3,Output_Sequence_of_Symbols_3] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(3));
[Output_Sequence_of_Bits_4,Input_Sequence_of_Symbols_4,Output_Sequence_of_Symbols_4] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(4));
[Output_Sequence_of_Bits_5,Input_Sequence_of_Symbols_5,Output_Sequence_of_Symbols_5] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(5));
[Output_Sequence_of_Bits_6,Input_Sequence_of_Symbols_6,Output_Sequence_of_Symbols_6] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(6));
[Output_Sequence_of_Bits_7,Input_Sequence_of_Symbols_7,Output_Sequence_of_Symbols_7] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(7));
[Output_Sequence_of_Bits_8,Input_Sequence_of_Symbols_8,Output_Sequence_of_Symbols_8] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(8));
[Output_Sequence_of_Bits_9,Input_Sequence_of_Symbols_9,Output_Sequence_of_Symbols_9] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(9));
[Output_Sequence_of_Bits_10,Input_Sequence_of_Symbols_10,Output_Sequence_of_Symbols_10] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(10));
[Output_Sequence_of_Bits_11,Input_Sequence_of_Symbols_11,Output_Sequence_of_Symbols_11] = PAM(Input_Sequence_of_Bits,M,Lb,Tc,Tsample,Tsymbol,SNR(11));

% Obtain the ber ( bit error rate ) for each SNR
ber1 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_1,Lb);
ber2 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_2,Lb);
ber3 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_3,Lb);
ber4 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_4,Lb);
ber5 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_5,Lb);
ber6 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_6,Lb);
ber7 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_7,Lb);
ber8 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_8,Lb);
ber9 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_9,Lb);
ber10 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_10,Lb);
ber11 = BER(Input_Sequence_of_Bits,Output_Sequence_of_Bits_11,Lb);

ber = [ber1 ber2 ber3 ber4 ber5 ber6 ber7 ber8 ber9 ber10 ber11];

% Obtain the ser ( symbol error rate ) for each SNR
ser1 = SER(Input_Sequence_of_Symbols_1,Output_Sequence_of_Symbols_1,Lb,M);
ser2 = SER(Input_Sequence_of_Symbols_2,Output_Sequence_of_Symbols_2,Lb,M);
ser3 = SER(Input_Sequence_of_Symbols_3,Output_Sequence_of_Symbols_3,Lb,M);
ser4 = SER(Input_Sequence_of_Symbols_4,Output_Sequence_of_Symbols_4,Lb,M);
ser5 = SER(Input_Sequence_of_Symbols_5,Output_Sequence_of_Symbols_5,Lb,M);
ser6 = SER(Input_Sequence_of_Symbols_6,Output_Sequence_of_Symbols_6,Lb,M);
ser7 = SER(Input_Sequence_of_Symbols_7,Output_Sequence_of_Symbols_7,Lb,M);
ser8 = SER(Input_Sequence_of_Symbols_8,Output_Sequence_of_Symbols_8,Lb,M);
ser9 = SER(Input_Sequence_of_Symbols_9,Output_Sequence_of_Symbols_9,Lb,M);
ser10 = SER(Input_Sequence_of_Symbols_10,Output_Sequence_of_Symbols_10,Lb,M);
ser11 = SER(Input_Sequence_of_Symbols_11,Output_Sequence_of_Symbols_11,Lb,M);

ser = [ser1 ser2 ser3 ser4 ser5 ser6 ser7 ser8 ser9 ser10 ser11];


%% Plot the BER diagram for SNR = 0:2:20 db
semilogy(SNR,ber); 

%% Plot the SER diagram for SNR = 0:2:20 db
semilogy(SNR,ser); 

