function identity = detector(Xf,num_frame_noise,threshold,mode)
%{

Xf: freq,time,channels of the stft of the full data.
num_frames_noise: number of frames (stft domain) that is considered noise
threshold: x times the noise floor (not in db,x times the energy)
mode = "binary" for 0 or 1 detection. "deviation" for energy deviation from
noise floor (energy- noisefloor)/noisefloor.  INCOMPLETE

%}

Nf_hat = Xf(:, 1:num_frame_noise, :);


identity = zeros(size(Xf(1,num_frame_noise+1:end,1)));
for i =1:size(Nf_hat,1)

  
  nf_k =  squeeze(Nf_hat(i, :, :));
  nf_k=nf_k';
  Rn_k = nf_k*nf_k'/size(nf_k, 2);

  xk_initial= squeeze(Xf(i,num_frame_noise+1:end,:))';
  Rx_k_initial= xk_initial*xk_initial'/size(xk_initial,2);

 [Vxn, Dxn, Wxn] = eig(Rx_k_initial, Rn_k);
 [eigs,idx] = sort(diag(Dxn),'descend');
  
 A_look = Vxn(:,idx(1));
 A_null=Vxn(:,idx(2:end));
    
 C = horzcat(A_look,A_null);
 f = vertcat(ones(1),zeros(1,size(C,2)-1)');

 w = pinv(Rx_k_initial)*C*pinv(C'*pinv(Rx_k_initial)*C)*f;
 
 noise_floor = norm(w'*nf_k)/length(nf_k);
 
 identity_k = zeros(1);
 for j =1:length(xk_initial)
   
    x_j = xk_initial(:,j);
    energy = norm(w'* x_j);
    
    if mode =="binary"
        if energy>= threshold *noise_floor
          identity_k = horzcat(identity_k,1);
        else
         identity_k = horzcat(identity_k,0);
        end
    elseif mode == "deviation"
        identity_k = horzcat(identity_k,(energy-noise_floor)/noise_floor);

    end
    
    

 end
identity_k=identity_k(1,2:end);
identity= vertcat(identity,identity_k);


end
identity = identity(2:end,:);




end