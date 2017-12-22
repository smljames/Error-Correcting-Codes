function Lambda = BerlekampMassey(S, t)
    % r: received signal
    % alpha: primitive element of GF(2^4)
    % t: error correcting capability

    % ----- initial condition -----
    r = 0;
    Lambda = [1 zeros(1, t)]; % Lambda = 1
    delta = 0;          % discrepency
    L = 0;
    B = [1 zeros(1,t)]; % B = 1
    
    % ----- Berlekamp-Massey Algorithm -----
    
    % --- find error position polynomial ---
    while(r ~= 2*t)
        r = r+1;
        delta = 0;
        for i = 0:L
            delta = delta + S(r-i)*Lambda(i+1);
        end
        if delta == 0
            B = [0, B(1:end-1)];
        else
            T = Lambda - delta.*[0 B(1:end-1)];
            if 2*L <= r-1
                B = Lambda./delta;
                Lambda = T;
                L = r-L;
            else
                Lambda = T;
                B = [0 B(1:end-1)];
            end
        end
    end
end