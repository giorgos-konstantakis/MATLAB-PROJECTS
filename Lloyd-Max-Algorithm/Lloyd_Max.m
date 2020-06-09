function [xq,centers,D,K_max,SQNR,entropy,c] = Lloyd_Max(x,N,x_min,x_max,epsilon)

% This function applies the Lloyd-Max algorithm in a signal , for example in
% a sound source or image source .
%
% The inputs of the function are :
% x : the input signal in vector form
% N : the number of bits used for modulation
% x_min : the minimum value for the normalized signal
% x_max : the maximum value for the normalized signal
% epsilon : the convergence condition for the algorithm
%
% The output of the function are :
% xq : the processed signal
% centers : the final centroids of the quantization
% D : the Distortion log
% K_max : the number of iterations of the algorithm
% SQNR : the SQNR value of the signal
% entropy : the entropy of the signal
% c : the probabilities of the quantizations's level values


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
    
   err_val = 0;
   c=zeros(length(x0));
   cond_mean=zeros(length(x0));
   
   % Calculate the quantization's zones
   Tk(1)=x_min;
   for i=2:(length(x0)-2)
       Tk(i) = (x0(1,i) + x0(1,i+1))/2 ;
   end
   Tk(i+1)=x_max;
   clear i
   
   % Calculate the quantizised signal xq
   % Also , calculating the the quantization's error and the conditional
   %  mean , so to find the  mean distortion due to quantization 
   for i=1:length(x(:,1))
       for z=1:(length(x0)-1)
           if (Tk(z)<x(i) && Tk(z+1)>=x(i))
               err_val=err_val+abs(x0(1,z+1)-x(i));
               cond_mean(z)=cond_mean(z)+x(i);
               c(z)=c(z)+1;
               break;
           end
       end
       xq(i,1)=x0(1,z); 
   end
   clear i
   clear z
   
   % Calculate the signal's distortion
   D(1,j) = err_val/length(x(:,1));
   
   % Calculate the new centroids
   for i=2:(length(x0)-2)
      if (c(i-1)~=0)
          x0(i)=cond_mean(i-1)/c(i-1);
      end
   end
   clear i
   
   % Calculating the SQNR
   S_Q_N_R = mean(x.^2)/mean((x-xq).^2) ; 
   SQNR(1,j) = 10 * log10(S_Q_N_R) ; % SQNR in dB
   
   % Check the distortion level , so to break or not the loop
   if( j>1 && abs(D(1,j)-D(1,j-1))<epsilon)
       K_max=j;
       break;
   end
   
end

% Calculation of Entropy in the output of the quantizer
c=c(1:(length(x0)-2),1);
entropy = 0;
for j = 1:(length(x0)-2)
    if(c(j)~=0)
        entropy = entropy - (c(j)/length(x(:,1)))*log2((c(j)/length(x(:,1)))) ;
    end
end
 c=c/length(x(:,1));
clear j
centers=x0(2:(length(x0)-1));
end