function [M, Rain] = CSC_ADMM_CPU( X, F, lambda, MaxIter)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
max_rho = 1e5; rho=1; tol = 1e-10;
[h,w] = size(X); K=size(F,3);
iter = 1; Cond = 1;

lambda=repmat(reshape(lambda,[1,1,K]),[h,w,1]);
for k=1:K
    Filters(:,:,:,k) = psf2otf(rot90(F(:,:,k),2),[h,w,n]);
end
FX = fft2(single(X));
MU = zeros(size(Filters));
Y = zeros(size(Filters));

C_Filters = conj(Filters);            % [h,w,K]
FTX = C_Filters.*repmat(FX,[1,1,K]);  % [h,w,K]
FTF  = sum(C_Filters.*Filters,3);     % [h,w]

while(iter < MaxIter && Cond)
    FR = FTX+rho*fft2(Y-MU);            % [h,w,K]
    FM = (FR-repmat(sum(Filters.*FR,3)./(rho+FTF), [1,1,K]) .*C_Filters)./rho;  % [h,w,K]    
    M = real(ifft2(FM));   % [h,w,K]
    Y = max(M+MU-lambda/rho,0);  % [h,w,K]
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
M = double(Y);
Rain = double(sum(real(ifft2(FM.*Filters)),3));    


