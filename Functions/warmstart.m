function model = warmstart(X_Fold,model)
% Use reweighted L2 matrix factorization to calculate U
 disp('Warm-start...');
 disp('initializing U...');
 [h,w,n]=size(X_Fold);
 X=reshape(X_Fold,[h*w,n]);
%% PCA and L1MF for warmstart
[~, U, oriV] = inexact_alm_rpca(X);
oriV = gather(oriV);
r = model.r;
U = gather(U(:,1:r));
V = oriV(:,1:r);
L = U*V';
%% Calculating the model
% noi = max(X-X(:, [2:n,1]), 0);
noi=X-L;
model.Sigma = var(noi(:))*50*0.05; %a simple initialization of mog parameters
model.D = InitDictNN( reshape(noi, [h,w,n]), model.f_size, model.video_path);     % [f_size,f_size,K]

%% Calculating A and B
m=size(X,1); 
A=zeros(r,r,m);B=zeros(r,m);
for i=1:m
    A(:,:,i)=(V'*V+0.001*eye(r))^-1*0.02;
    B(:,i)=V'* (L(i,:))'/0.02;
end
disp('Warm start is over. ');
%% output
model.A=A; model.B=B;
model.U=U;
model.v=repmat(mean(oriV(:,1:r)),[model.batch_size,1]);
model.r=r;
%% update dictionry online
model.FA=zeros(length(model.f_size));
model.FB=zeros(max(model.f_size)^2,length(model.f_size));
%% our update dictionry online 
% model.FA=zeros(length(model.f_size)*max(model.f_size)^2);
% model.FB=zeros(length(model.f_size)*max(model.f_size)^2,1);
end

function  Filters  = InitDictNN( R, f_size, video_path )
% pca for initialize D
FInd = FilterInd(f_size);
maxf_size = max(f_size);
K = length(f_size);
Patches = Video2Patch( R, maxf_size);
Patches = Patches- mean(Patches,2); 
temp = Patches*Patches';
[U, ~, ~] = svd(temp);
D = Normalize(max(U(:,2:K+1),0)); 
Filters = zeros(size(D)); Filters(FInd) = D(FInd);
display_filters=1;
if display_filters==1
    Filters_Fold = reshape(Filters, [maxf_size, maxf_size*K]);
    Filters_Fold = imresize(Filters_Fold, 10);
    imshow(Filters_Fold*10)
    imwrite(Filters_Fold*10, [video_path, 'filters.jpg'])
end
end
