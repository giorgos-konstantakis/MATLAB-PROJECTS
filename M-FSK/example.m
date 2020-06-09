clc ; clear all ;
% % % Example of M-FSK algorithm % % %

% Number of bits
Lb = 10000;

% Input sequence
Input = randsrc(Lb,1,[0 1]);

% M_PSK's order
M=2; % choose  M=4(one symbol per 2 bits) or M=8(one symbol per 3 bits) or M=16(1 symbol per 4 bits)

% Time period of symbol
Tsymbol=40; 
% Time period of cosine's carrier
Tc=4;
% Time period of sampling
Tsample=1;

% Coding bits
% Set 0 for normal code or 1 for gray code
coding = 0;


% SNR in db ( we set many SNR values , to calculate the BER diagram )
SNR = 0;

% Energy of each symbol ( we usually normalize it to 1 )
Es = 1;

% M-PSK algorithm execution
[output_sequence_of_bits,input_sequence_of_bits,q]=M_FSK(Input,M,Lb,Tc,Tsample,Tsymbol,SNR,Es,coding);

% you can check the error with biterr
error = biterr(output_sequence_of_bits,input_sequence_of_bits);
error = error/length(output_sequence_of_bits(:,1))
