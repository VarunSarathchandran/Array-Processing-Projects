function Xf_filtered = updated_filtering(DetectedLocations, Xf, Nf_hat, alpha)

factor = 1;
[num_freqbd, ~] = size(DetectedLocations);
num_frame = size(Xf, 2);
Sf_hat = zeros(num_freqbd, num_frame);

for i = 1:num_freqbd
    LocationsPerFreq = DetectedLocations(i, :);
    XfPerFreq = squeeze(Xf(i, :, :)).';
    nf = squeeze(Nf_hat(i, :, :));
    nf = nf.';
    Rn = nf*nf'/size(nf, 2);
    Rx = Rn;
    for j = 1:num_frame
        Rx = Rx*alpha + XfPerFreq(:, j)*XfPerFreq(:, j)'*(1-alpha);
        if LocationsPerFreq(j) == 0
            Rn = Rn*alpha + XfPerFreq(:, j)*XfPerFreq(:, j)'*(1-alpha);
        else
            Rn = Rn;
        end
        [Vxn, Dxn] = eig(Rx, Rn);
        [eigs,idx] = sort(diag(Dxn),'descend');
        Uxn = Vxn(:,idx);
        Qxn = inv(Uxn');
        U1 = Uxn(:, 1);
        Q1 = Qxn(:, 1);
        eig_s = eigs(1)-1;
        %w = eig_s*(U1)/(eig_s + factor);
        % w = eig_s*a/(eig_s + factor);
        w = (eig_s/(eig_s+factor))*pinv(Rn)*Q1;
        Sf_hat(i, j) = w'*XfPerFreq(:, j);
    end
end
Xf_filtered = Sf_hat;



