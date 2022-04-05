function [K,DVW] = DWVT(x)

%Discrete Wigner-Ville Distribution
% X = 1D signal amplitude
% K = Kernels of WVD
% DVW = Discrete Wigner-Ville Distribution


if (rem(length(x),2) == 0)
    error("Window Size must be in odd number")
end

N = length(x)-1;
X=fft(x);
X=[X(1:N/2+1);zeros(N,1);X(N/2+2:N+1)];
x=2*ifft(X);
x=[zeros(N,1);x;zeros(N,1)];
z = hilbert(x);
X=zeros(2*N+1);

for k=1:2*N+1  
    X(:,k)=z(k+(0:2*N)); 
end
 K=X.*conj(flipud(X));
 DVW=fft(K([N+1:2*N+1,1:N],:)); 
end