function [M, Rain, Rains] = CSC_ADMM_GPU( X, F, lambda, MaxIter)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
max_rho = 1e5; rho=1; tol = 1e-10;
[h, w, n_frames] = size(X); K=size(F,3);
iter = 1; Cond = 1;

lambda=gpuArray(repmat(reshape(lambda,[1,1,1,K]),[h, w,n_frames, 1]));
F = gpuArray(single(F));
Filters = gpuArray.zeros(h,w,n_frames,K);
for k=1:K
    Filters(:,:,:,k) = psf2otf(rot90(F(:,:,k),2),[h,w,n_frames]);
end
FX = gpuArray(fft2(single(X)));
MU =  gpuArray.zeros(size(Filters),'single');
Y =  gpuArray.zeros(size(Filters),'single');

C_Filters = conj(Filters);            % [h,w,n_frames, K]
FTX = C_Filters.*repmat(FX,[1,1,1,K]);  % [h,w,n_frames,K]
FTF  = sum(C_Filters.*Filters,4);     % [h,w,n_frames]

while(iter < MaxIter && Cond)
    FR = FTX+rho*fft2(Y-MU);            % [h,w,n_frames,K]
    FM = (FR-repmat(sum(Filters.*FR,4)./(rho+FTF), [1,1,1,K]) .*C_Filters)./rho;  % [h,w,n_frames,K]   
    M = real(ifft2(FM));   % [h,w,n_frames,K]
    Y = max(M+MU-lambda/rho,0);  % [h,w,n_frames,K]
    MU = MU + M - Y;
    if(rho < max_rho)
        rho = rho*1.01;
    end
    if(mod(iter,10)==0)
        ConvergenceError = mean((M(:)-Y(:)).^2);
        Cond = (ConvergenceError>tol);
    end
    iter = iter+1;
end
M = double(gather(Y));
Rains = double(gather(real(ifft2(FM.*Filters))));
% Rain = sum(Rains,4);
Rain = max(sum(Rains,4)-5e-3, 0);
clear F Filters lambda FX MU Y




