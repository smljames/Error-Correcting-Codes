function err_location = ErrorLocation(alpha, Lambda, m, t)
    k = 1;
    err_location = gf(zeros(1, t),m); % initial the error position vector
    
    % --- exhaustive search to find the root of polynomial ---
    for i = 0:2^m-2 % alpha^0 --> alpha^(2^(m-2))
        tmp = 0;
        for j = 1:length(Lambda)
            tmp = tmp + Lambda(j)*(alpha^i)^(j-1);
        end
        if tmp == 0
            err_location(k) = 1/(alpha^i);
            k = k+1;
        end
    end
end