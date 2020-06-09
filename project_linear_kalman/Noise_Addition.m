function [w,v] = Noise_Addition(Qk,Rk,n)

% This function adds process noise and observation noise to the plant

% Inputs
% Qk = process noise variance
% Rk = observation noise variance
% n = time vector's length

% w = process noise ( Qk is process's noise variance )
w = sqrt(Qk)*randn(n,1);

% v = observation noise ( Rk is observation's noise variance )
v = sqrt(Rk)*randn(n,1);

end