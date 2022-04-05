function rc = arburg(x,L)

% rc =  reflection coefficient
% x = segment/kernels values
% L = Order

if ~iscell(x)
    if size(x,1) == 1, x = x'; 
    end
    x = {x};
end
rc(1) = 1;
ns = length(x);
f = x; 
b = x;

for j = 1:L
    num = 0; den = 0;
    v1 = cell(ns,1); v2 = cell(ns,1);
    for i = 1:ns
        v1{i} = b{i}(1:size(b{i},1)-1,:); 
        v2{i} = f{i}(2:end,:);
        num = num+sum(dot(v1{i},v2{i}));
        den = den+sum(dot(v1{i},v1{i}))+sum(dot(v2{i},v2{i}));
    end  
    k = -2*num/den;
    rc = [rc k];
    for i = 1:ns
        f{i} = v2{i}+k*v1{i}; b{i} = v1{i}+k*v2{i};
        %var1(i) = 1/2*((f{i}).^2 +(b{i}).^2 );
    end %for i = 1:ns,
end %for j = 1:L,
end
