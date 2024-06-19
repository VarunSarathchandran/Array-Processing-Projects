function x = gendata_conv(s, P, N, sigma)

x = zeros(1, N*P);

for t = 1:N*P
    for k = 1:N
        t_k = (t-1)/P - (k-1);
        if (t_k>=0 && t_k<0.25)||(t_k>=0.5 && t_k<0.75)
            h = 1;
        elseif (t_k>=0.25 && t_k<0.5)||(t_k>=0.75 && t_k<1)
            h = -1;
        else
            h = 0;
        end
        x(t) = x(t) + h*s(k);
    end
end
n_real = randn(1, N*P);
n_imag = randn(1, N*P);
n_real = sigma/sqrt(2) * n_real;
n_imag = sigma/sqrt(2) * n_imag;
n = n_real + 1i*n_imag;
x = x + n;

end

