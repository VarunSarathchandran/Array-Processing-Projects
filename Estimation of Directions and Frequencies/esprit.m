function theta_hat = esprit(X, d)

[M, ~] = size(X);
Z_upper = X(1:M-1, :);
Z_lower = X(2:M, :);
[Uz, ~, ~] = svd([Z_upper; Z_lower]);

ux = Uz(1:(M - 1), 1:d);
uy = Uz(M:end, 1:d);

if (M - 1) < d
    phi_full = ux' * inv(ux * ux') * uy;
else
    phi_full = inv(ux' * ux) * ux' * uy;
end

[~, phi] = eig(phi_full);

theta_hat = zeros(d,1);

for i = 1:d
    phase = angle(phi(i, i));
    %if phase < 0
    %    phase = phase + 2*pi;
    %end
    theta_hat(i) = asind(phase/(2*pi*0.5));
end

theta_hat = sort(theta_hat);

end