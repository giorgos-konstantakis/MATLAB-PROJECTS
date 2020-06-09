function xq = ADM(y,K,Delta,fs,fd)

% This function apply Adaptive Delta Modulation ( ADM ) in a signal , for
% example a sound source or an image source . 
% The inputs of the function are :
% y : the signal to be processed in vector form
% K : the changing value of step Delta ( good choice is = 1.5 )
% Delta : the initial Delta step ( choose a small value , around zero )
% fs : the new sampling frequency ( fs>fd and fs>>Nyquist_Frequency , in order to have enough data to perform the ADM )
% fd : the original sampling frequency
% The output of the function is :
% xq : the processed signal 

M = round(fs/fd);

% Oversampling the signal
x = interp(y,M);

% Initialization of algorithm's matrices
N = length(x(:,1));
e = zeros(N+1,1);
delta = zeros(N+1,1);
b = zeros(N+1,1);
eq = zeros(N+1,1);
xq = zeros(N+1,1);
xq_2 = zeros(N+1,1);
eq_2 = zeros(N+1,1);
delta(1,1)=Delta;

% ADM algorithm
for i=2:(N)
    
    % % % MODULATOR % % %
    
    % Quantizer
    if abs((x(i,1)-xq(i-1,1))+1)<=abs((x(i,1)-xq(i-1,1))-1)
        b(i,1)=-1;
    else
        b(i,1)=1;
    end
    
    % Logic decision
    if b(i,1) == b(i-1,1)
        delta(i,1) = delta(i-1,1)*K;
    else
        delta(i,1) = delta(i-1,1)/K;
    end
    
    % Multiplier
    eq(i,1) = delta(i,1)*b(i,1);
    
    % Summer
    xq(i,1) = eq(i,1) + x(i-1,1);
    
    % % % DEMODULATOR % % %
    
    % Multiplier
    eq_2(i,1) = b(i,1)*delta(i,1);
    
    % Summer
    xq_2(i,1) = eq_2(i,1) + xq_2(i-1,1);

end

% Output signal
xq = xq_2(2:(N+1));

end