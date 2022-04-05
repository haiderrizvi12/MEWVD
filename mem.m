function [K1,DVW] = mem(K,size,len,order,fs,clength,sig)

%Pre-release version v1.0.1
% Modified Kernels
% Modied DVW 
% K = Kernels values 
% size = Window size
% len = Maximum modified kernels
% order = AR order
% fs = sampling frequency
% clength = contour plot maximum length (optional)
% sig  = original signal (Optional)


if (rem(size,2) == 0)
    error("Window Size must be in odd number")
end

K1 = zeros(length(K));

L= length(K(:,1));

Ns = (L-1)/2;
J = (size-1)/2;
rs = Ns+1-J;
re = Ns+1+J;
req = Ns- J;
loop = 0;
for i=1:len:L
 if((loop+len)<L) 
for j =1:1:len
    %x = K(rs:1:re,loop+j);
    x(:,j) = K(rs:1:re,loop+j);  
    %a = arburg(x,order);
    %K1(:,j) = Proj(x,a,req,'post'); 
end
%a = arburg(x,order);
ref = arburg_my(x,order);
a = rc2lpc(ref(2:end));
% c = mean(x,2);

 % a=arburg(c,order);
K1(:,i:loop+len) = Proj(x,a,req,'post'); 
 loop = loop+len;
 end
end

if (length(K1(1,:)) < length(K(1,:)))
    
    x = K(rs:1:re, length(K1(1,:))+1:1:L);
    ref = arburg_my(x,order);
a = rc2lpc(ref(2:end));
 
    K1(:,length(K1(1,:))+1:1:L) = Proj(x,a,req,'post'); 
end
    
K1(rs:Ns+1,1:L) = K(rs:Ns+1,1:L);


%for i=Ns+1-J+1:1:Ns+1
%K1(i,:) = K(i,:);
%end
h = K1(1:Ns,:);
 h = conj(flip(h));
 %u=0;
%for i = Ns+2:1:2*Ns
 %   u = u+1;
    K1(Ns+2:2*Ns+1,:) = h;
   K1 = K1.* g;
%endca
DVW=fft(K1([Ns+2:2*Ns,1:Ns+1],:));
DVW = real(DVW);
%DVW=real(DVW([Ns+2:2*Ns+1,1:Ns+1],:));
%if nargout < 1

figure
if (nargin >= 7)
if (rem(length(sig),2) == 0)
    error("Input Signal must be in odd number")
end
I = DVW(1:clength,:);
t = 1:length(I);

N = length(sig)-1;

X=fft(sig);
X=[X(1:N/2+1);zeros(N,1);X(N/2+2:N+1)];
x=2*ifft(X);



tfrqview(real(I),x,t,'tfrspwv');
else
    N = (length(DVW)-1)/2;
t = (0:N)/fs;
 f = (0:N)'/(N+1) * fs;
f = linspace(0,f(end),N);
imagesc(t*1e6,f/1e6,real(DVW));
set(gca,'YDir','Normal')

l = caxis;
caxis([0 l(2)]);

end
title(size);

end

%end


% L= length(K(:,1));
% 
% Ns = (L-1)/2;
% J = (size-1)/2;
% for i=1:1:Ns
% 
%     x = K(i,Ns-J:1:Ns+J);
%     x=x';
%  K1(:,i)= lpredict(x,order,Ns-J,'pre');
% 
%  K2(:,i) = lpredict(x,order,Ns-J,'post');
% 
% 
% end
%  K1 = K1';
%   K2 = K2';
% for i=Ns-J+1:1:Ns+J+1
% K1(:,i) = K(1:Ns,i);
% end
% K1(:,Ns+J+2:1:L) = K2;
% 
% h = K1(1:Ns,:);
%  h = conj(flip(h));
%  %u=0;
% %for i = Ns+2:1:2*Ns
%  %   u = u+1;
%     K1(Ns+2:2*Ns+1,:) = h;
% %end
% DVW=fft(K1([Ns+2:2*Ns,1:Ns+1],:));
% figure
% imagesc(real(DVW));
% set(gca,'YDir','Normal')


