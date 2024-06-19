function [X, A, S] = gendata(M, N, Delta, theta, f, SNR)
%GANDATA Summary of this function goes here
%   Detailed explanation goes here
%   Currently we set the antenna response to 1

d = length(f);

A = zeros(M, d);
S = zeros(d, N);

theta = theta*pi/180;

for x = 1:M
    for y = 1:d
        A(x, y) = exp(2i*pi*(x - 1)*Delta*sin(theta(y)));
    end
end

for x = 1:d
    for y = 1:N
        S(x, y) = exp(2i*pi*f(x)*(y - 1));
    end
end


if SNR == inf
    X = A*S;
else
    mu = 0;
    sigma = sqrt(power(10, -0.1*SNR))/sqrt(2);

    N_R = mu + sigma * randn(M,  N);
    N_C = mu + sigma * randn(M,  N);
    N = N_R + N_C*1i;

    X = A*S + N;
end

