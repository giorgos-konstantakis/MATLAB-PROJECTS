function ber = BER(x_input,x_output,Lb)

% This function calculates the bit error rate ( BER ) of an M-PAM
% modulation.
% 
% The function's inputs are :
% x_input : the input sequence of bits to be transmitted
% x_output : the output sequence of bits
% Lb : the number of the transmitting bits
%
% The function's outputs are :
% ber : the bit error rate ( SER )

ber=0;
for i=1:Lb
    if x_input(i,1)~=x_output(i,1)
        ber = ber +1;
    end
end
    
ber = ber/Lb;

end
