
function A = rc2lpc(K)
% RC2LPC Convert Reflection Cofficients to Linear Prediction Coefficients.


        n = 1; 
        for i = 2:length(K) %order
            a = n;
            for j = 1:(i-1)
                n(j) = a(j) + K(i)*a(i-j);
            end
            n(i) = K(i);
        end
        A = [1, n];
    
end
