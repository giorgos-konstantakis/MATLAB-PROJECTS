clc; clear all;

% Basic System
% Here we create our system , in which the Kalman Filter is going to be
% implemented

% Choose G as you wish
G=tf( [ 1 100000 100 ] , [ 1 128 4481 30540 75900 70000 ] );
sys1=ss(G);
sys=c2d(sys1,1);
A=sys.A;
B=sys.B; 
C=sys.C;
e=eig(G);
t=(0:100);
n=length(t);
u=ones(n,1);

% Covariance of process and observation noise
% Choose the values of Qk , Rk 
Qk=1;
Rk=1;

% Process and Observation noise
[w,v] = Noise_Addition(Qk,Rk,n);

% Model's output
y = lsim(sys,u);

% System's output with process noise
y_sys=lsim(sys,u+w);

% System's output with process noise and observation noise
y_meas=y_sys+v;

% Linear Kalman Filter's implementation to our system
[X_new,y_filter,Kgain,errcov,P] = Linear_Kalman_Filter(Qk,Rk,n,A,C,y_meas);


