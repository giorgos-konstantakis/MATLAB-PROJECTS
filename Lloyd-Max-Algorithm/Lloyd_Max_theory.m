function [xq,centers,D,K_max,c] = Lloyd_Max_theory(x,N,x_min,x_max,epsilon,s2,mean)

% This function applies the Lloyd-Ma algorithm in a signal , for example in
% a sound source or image source .
%
% The inputs of the function are :
% x : the input signal in vector form
% N : the number of bits used for modulation
% x_min : the minimum value for the normalized signal
% x_max : the maximum value for the normalized signal
% epsilon : the convergence condition for the algorithm
% s2 : the variance value of the signal's pdf
% mean : the mean value of signal's pdf
%
% The output of the function are :
% xq : the processed signal
% centers : the final centroids of the quantization
% D : the Distortion log
% K_max : the number of iterations of the algorithm
% c : the probabilities of the signal's values

syms xx
syms xxx

% number of quantization's levels
L = 2^N;

% initial quantization's step ( for omogenous quantization )
h = (x_max-x_min)/L;

% Creating the initial values of quantization's levels
c1 = -L/2+1;
c2 = L/2-(L/2-1);
x0=zeros(1,L);
for i = 1:L
    if(i <= L/2)
        x0(1,i) = (c1-1/2) * h;
        c1=c1+1;
    else
        x0(1,i) = (c2-1/2) * h;
        c2=c2+1;
    end
end
clear i
clear c1
clear c2
clear h

x0=[x_min x0 x_max]; 
xq=x;
% Loop of the algorithm
for j=1:10000

   c=zeros(length(x0));
   
   % Calculate the quantization's zones
   Tk(1)=x_min;
   for i=2:(length(x0)-2)
       Tk(i) = (x0(1,i) + x0(1,i+1))/2 ;
   end
   Tk(i+1)=x_max;
   clear i
   
   % Calculate the quantizised signal xq 
   for i=1:length(x(:,1))
       for z=1:(length(x0)-1)
           if (Tk(z)<x(i) && Tk(z+1)>=x(i))
               c(z)=c(z)+1;
               break;
           end
       end
       xq(i,1)=x0(1,z); 
   end
   clear i
   clear z
   
   % Calculate the signal's distortion with normal distribution f(x)
   D(1,j)=0;
   for i=2:length(x0)-1
       fx(xx) = (1/(sqrt(2*pi*s2)))*exp(-((xx-mean)^2)/(2*s2));
       expr(xx) = ((xx - x0(1,i))^2)*fx(xx);
       a = Tk(i);
       b = Tk(i-1);
       F = vpaintegral(expr,b,a);
       F = double(F);
       D(1,j) = D(1,j)+F;
   end
   clear i
   
   % Calculate the new centroids
   for i=2:length(x0)-1
       fx(xxx) = (1/(sqrt(2*pi*s2)))*exp(-((xxx-mean)^2)/(2*s2));
       expr_1(xxx) = xxx*fx(xxx);
       I1 = vpaintegral(expr_1,xxx,[Tk(i-1) Tk(i)] );
       I2 = vpaintegral(fx,xxx,[Tk(i-1) Tk(i)] );
       I1 = double(I1);
       I2 = double(I2);
       x0(1,i)=I1/I2;
   end
   clear i
   
   % Calculating the SQNR
   %S_Q_N_R = mean(x.^2)/mean((x-xq).^2) ; 
   %SQNR(1,j) = 10 * log10(S_Q_N_R) ; % SQNR in dB
   
   % Check the distortion level , so to break or not the loop
   if( j>1 && abs(D(1,j)-D(1,j-1))<epsilon)
       K_max=j;
       break;
   end
   
end

% Calculation of Entropy in the output of the quantizer
c=c(1:(length(x0)-2),1);
c=c/length(x(:,1));

clear j
centers=x0(2:(length(x0)-1));
end