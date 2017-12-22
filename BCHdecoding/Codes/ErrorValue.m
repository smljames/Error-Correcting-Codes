function err_value = ErrorValue(Z, Lambda, err_location, m, t)
    n = length(err_location);
    err_value = gf(zeros(1, t), m);
    for i = 1: n % err_value(1) --> err_value(n)
        if err_location(i)~=0
            tmp1 = 0;
            tmp2 = Lambda(1) * err_location(i);
            for j = 1:length(Z)
                tmp1 = tmp1 + Z(j)*err_location(i)^(-1)^(j-1);
            end
            for j = 1:n
                if j ~= i
                    tmp2 = tmp2 * (1 - err_location(j)*err_location(i)^(-1));
                end
            end
            err_value(i) = -tmp1/tmp2;
        end
    end