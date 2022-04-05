function [y, a]=Proj(x, a, npred, pos)

if nargin<3
    error('Not enough input arguments');
end

if nargin<4
    pos='post';
end
cols=size(x,2);
if cols>1
    y=zeros(npred, cols);
    
    for k=1:size(x,2)
        [y(:,k) ]=Proj(x(:,k), a, npred, pos);
    end
    return

if nargin==4 && strcmpi(pos,'pre')
    x=x(end:-1:1);
end


cc=-a(2:end);

y=zeros(npred,1);
y(1)=cc*x(end:-1:end-np+1);
for k=2:min(np,npred)
    y(k)=cc*[y(k-1:-1:1); x(end:-1:end-np+k)];
end
for k=np+1:npred
    y(k)=cc*y(k-1:-1:k-np);
end
if nargin==4 && strcmpi(pos,'pre')
    y=y(end:-1:1);
end

return
end