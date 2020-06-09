function ser = SER(x_input,x_output,Lb,M)

% This function calculates the symbol error rate ( SER ) of an M-PAM
% modulation.
% 
% The function's inputs are :
% x_input : the input sequence of bits to be transmitted
% x_output : the output sequence of bits
% Lb : the number of the transmitting bits
% M : the order of the M-PAM ( M = 4 or M = 8 )
%
% The function's outputs are :
% ser : the symbol error rate ( SER )

ser=0;
for i=1:(Lb/log2(M))
    if x_input(i,1)~=x_output(i,1)
        ser = ser +1;
    end
end
    
ser = ser/(Lb/log2(M));

end
