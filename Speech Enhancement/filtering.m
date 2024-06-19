function Xf_filtered = filtering(DetectedLocations, Xf, FrameLength, num_frames, num_subframes, Nf_hat)

factor = 1;
[num_freqbd, ~] = size(DetectedLocations);

Sf_hat = zeros(num_freqbd, num_subframes);
for i = 1:num_frames
    Locations = DetectedLocations(:, (i-1)*FrameLength+1:i*FrameLength);
    Xf_tmp = Xf(:, (i-1)*FrameLength+1:i*FrameLength, :);

    for j = 1:num_freqbd
        LocationsPerFreq = Locations(j, :);
        xf_tmp = squeeze(Xf_tmp(j, :, :));
        xf_tmp = xf_tmp.';
        idx_Rx = find(LocationsPerFreq == 1);
        idx_Rn = find(LocationsPerFreq == 0);
        if isempty(idx_Rx)
            continue;
        else 
            % if (isempty(idx_Rn))||(length(idx_Rn)<4)
            %     nf = squeeze(Nf_hat(j, :, :));
            %     nf = nf.';
            %     Rn = nf*nf'/size(nf, 2);
            % else
            %     Rn = xf_tmp(:, idx_Rn)*xf_tmp(:, idx_Rn)'/size(xf_tmp(:, idx_Rn), 2);
            % end
            nf = squeeze(Nf_hat(j, :, :));
            nf = nf.';
            Rn = nf*nf'/size(nf, 2);
            Rx = xf_tmp(:, idx_Rx)*xf_tmp(:, idx_Rx)'/size(xf_tmp(:, idx_Rx), 2);
            [Vxn, Dxn] = eig(Rx, Rn);
            [eigs,idx] = sort(diag(Dxn),'descend');
            Uxn = Vxn(:,idx);
            Qxn = inv(Uxn');
            U1 = Uxn(:, 1);
            Q1 = Qxn(:, 1);
            eig_s = eigs(1)-1;
            w = eig_s*(U1)/(eig_s + factor);
            %w = (eig_s/(eig_s+factor))*pinv(Rn)*U1;
            % w = eig_s*a/(eig_s + factor);
            Sf_hat(j, (i-1)*FrameLength+1:i*FrameLength) = w'*xf_tmp;
        end
    end
end
Xf_filtered = Sf_hat;

