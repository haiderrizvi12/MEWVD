function A = rc2lpc(K)

% RC2LPC Convert Reflection Cofficients to Linear Prediction Coefficients.


p = length(K); % Order of linear prediction
        a1 = K(1); % a1(1) = k1
        for i = 2:p
            a = a1;
            for j = 1:(i-1)
                a1(j) = a(j) + K(i)*a(i-j);
            end
            a1(i) = K(i);
        end
        A = [1, a1];
    
end