function model = warmstart_trans(X_Fold,model)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
disp('Warm-start...');
 %% prealignment
disp('preAlignment...');
[height,wide,n] = size(X_Fold);
TX = X_Fold; Tau = zeros(model.type,n);
% [TX,Tau] = preAlign(TX,median(TX,3),Tau);
[TX,~] = preAlign(TX,TX(:,:,ceil(n/2)),Tau);

disp('initializing U...');
X = reshape(TX,[height*wide,n]);
%% PCA and L1MF for warmstart
[~, U, V] = inexact_alm_rpca(X);
r = min(model.r,size(U,2));
U = gather(U(:,1:r));
V = gather(V(:,1:r));
model.B_Fold=reshape(U*V(1,:)',[height,wide]); model.r=r;

%% Calculating the model
noi=X-U*V';
model.Sigma =var(noi(:))*50*0.05; %a simple initialization of gaussian parameters
model.D = InitDictNN(reshape(noi, [height,wide,n]), model);     % [f_size,f_size,K]

%% update filter iterative parameters
model.FA=zeros(length(model.f_size));
model.FB=zeros(max(model.f_size)^2,length(model.f_size));
disp('Warm start is over. ');
end

function  Filters  = InitDictNN( R, model )
% pca 
f_size=model.f_size;
FInd = FilterInd(f_size);
maxf_size = max(f_size);
K = length(f_size);
Patches = Video2Patch( R, maxf_size, model.use_gpu );
Patches = Patches- mean(Patches,2); 
temp = Patches*Patches';
clear Patches
[U, ~, ~] = svd(temp);
D = Normalize(max(U(:,2:K+1),0)); 
Filters = zeros(size(D)); Filters(FInd) = D(FInd);
end