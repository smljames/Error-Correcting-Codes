function err_value = Forney(S, Lambda, err_location, m, t)
    % S: syndrome vector
    % Lambda: error location polynomial
    % err_location: error location vector
    % t: error correcting capability
    n = length(err_location);
    Omega = conv(S, Lambda);
    Omega = Omega(1:2*t); % mod x^2t
    err_value = gf(zeros(1,n),m);
    
    for i = 1: n % err_value(1) --> err_value(n)
        if err_location(i)~=0
            tmp1 = 0;
            tmp2 = err_location(i);
            for j = 1:2*t
                tmp1 = tmp1 + Omega(j)*err_location(i)^(-1)^(j-1);
            end
            for j = 1:n
                if j ~= i
                    tmp2 = tmp2 * (1 - err_location(j)*err_location(i)^(-1));
                end
            end
            err_value(i) = tmp1/tmp2;
        end
    end
end