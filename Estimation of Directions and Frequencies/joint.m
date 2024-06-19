function [theta_hat,f_hat] = joint(X, d, m)

[M, N] = size(X);
jthresh = 1.0e-8;
N_col = N-m+1;

Z = zeros(m*M, N_col);

for i = 1:m
    Z((i-1)*M+1:i*M, :) = X(:, i:N_col+i-1);
end

[U_full, ~, ~] = svd(Z);
U = U_full(:, 1:d);

% Formulate T^-1 * phi * T, phi contains frequency information
U_phi_x = U(1:M*(m-1), :);
U_phi_y = U(M+1:end, :);

A_phi = pinverse(U_phi_x, d)*U_phi_y;

% Formulate T^-1 * theta * T, theta contains angle information
U_theta_x = zeros(m*(M-1), d);
U_theta_y = zeros(m*(M-1), d);

for i = 1:m
    U_theta_x((i-1)*(M-1)+1:i*(M-1), :) = U((i-1)*M+1:i*M-1, :);
    U_theta_y((i-1)*(M-1)+1:i*(M-1), :) = U((i-1)*M+2:i*M, :);
end

A_theta = pinverse(U_theta_x, d)*U_theta_y;

% Joint diagonalization
[~, D] = joint_diag([A_phi, A_theta], jthresh);

f_hat = zeros(d, 1);
theta_hat = zeros(d, 1);

for i = 1:d
    phase = angle(D(i, i));
    if phase < 0
        phase = phase + 2*pi;
    end
    f_hat(i) = phase/(2*pi);

    phase = angle(D(i, i+2));
    theta_hat(i) = asind(phase/(2*pi*0.5));
end

[theta_hat, sorting_indices] = sort(theta_hat);
f_hat = f_hat(sorting_indices);

end



