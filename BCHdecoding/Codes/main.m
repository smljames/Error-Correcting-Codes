% This code is to implement BCH decoding algorithm
% I use two algorithm:
%   1. Berlekamp-Massey Algorithm
%   2. Euclidean Algorithm

clc; clear;

% ----- set up Galois Field primitive element and field -----
t = 5;         % error correcting capability
m = 6;
alpha = gf(2,m); % primitive element of GF(2^m)

GF = gf(zeros(1,2^m),m);
GF(1) = 0;
for i = 0:2^m-2 % alpha^0 --> alpha^14
    GF(i+2)=alpha^i;
end


% ----- received signal -----
r = [0,0,alpha^6,0,0,0,0,alpha^20,alpha^3,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,alpha^62,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,alpha^15];




% ----- calculate syndrome -----
S = gf(zeros(1,2*t),m);          % initial syndrome vector
beta = gf(zeros(1,length(r)),m); % initial beta vector

for i = 1:length(r)
    beta(i) = alpha^(i-1); % beta = [a^0, a^1, ..., a^7]
end

for i = 1:2*t
    beta_ = beta.^i;       % beta_ = [(a^0)^i, (a^1)^i, ..., (a^7)^i]
    for j = 1:length(r)
        S(i) = S(i) + r(j)*beta_(j); % inner product(r, beta_)
    end
end



% ---------- part 1 ---------- %
%{
% output vector power: low to high
% Example:
%   err_location_alpha = 
%       51  35  8   7   2
%   It means that error location polynomial is
%   (1+a^51X)(1+a^35X)(1+a^8X)(1+a^7X)(1+a^2X)

% ----- Berlekamp-Massey Algorithm -----
Lambda = BerlekampMassey(S, t);
err_location = ErrorLocation(alpha, Lambda, m, t);
Lambda_alpha = alphapower(Lambda, GF);
err_location_alpha = alphapower(err_location, GF);

% ----- Forney Algorithm -----
err_value = Forney(S, Lambda, err_location, m, t);
err_value_alpha = alphapower(err_value, GF);
%}



% ---------- part 2 ---------- %
%{
% ----- Euclidean Algorithm -----

% --- initial conditions ---
a = gf([1, zeros(1, 2*t)], m);
S = fliplr(S);
f = gf([zeros(1, t) 0], m);
g = gf([zeros(1, t) 1], m);

[Z, Lambda] = Euclidean(a, S, f, g);

% --- reverse vectors ---
S = fliplr(S); % power: low to high
Z = fliplr(Z); % power: low to high
Lambda = fliplr(Lambda); % power: low to high

err_location = ErrorLocation(alpha, Lambda, m, t);
err_value = ErrorValue(Z, Lambda, err_location, m, t);
err_location_alpha = alphapower(err_location, GF);
err_value_alpha = alphapower(err_value, GF);
%}
