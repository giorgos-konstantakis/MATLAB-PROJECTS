function [X_new,y_filter,Kgain,errcov,P] = Linear_Kalman_Filter(Qk,Rk,n,A,C,y_meas)

% This function implements the Linear Kalman Filter in the system

% Inputs
% Qk = process noise variance
% Rk = observation noise variance
% n = time vector's length
% A = system's matrix A
% C = system's matrix C
% y_meas = real system's output (system with process and observation noise)

% Outputs
% X_new = the new states of the system after filter's processing
% y_filter = filter's output
% Kgain = Kalman Gain matrix
% errcov = filter's error
% p = error covariance matrix


% Initial error covariance
P = Qk;  

% Initial states' values
x = zeros(5,1);  

% Filter's output
y_filter = zeros(n,1);

% Filter's error 
errcov = zeros(n,1);

% New state values matrix
X_new=zeros(5,n);

for i = 1:n
    
  X_new(:,i)=x;
  
  % Measurement update
  
  % Kalman Filter's Gain = Kgain
  Kgain = P*C'/(C*P*C'+Rk);
  
  % New state's calculation
  x = x + Kgain*(y_meas(i)-C*x);   % x[n|n]
  
  % Here P matrix is calculated with Joseph's formula ( more robust )
  P = (eye(5)-Kgain*C)*P*(eye(5)-Kgain*C)'+Kgain*Rk*Kgain';      % P[n|n]
    
  % Filter's output
  y_filter(i) = C*x;
  
  % Calculation of filter's error
  errcov(i) = C*P*C';

  % Time update
  x = A*x ;        % x[n+1|n]
  P = A*P*A' + Qk;     % P[n+1|n]
  
end





















end