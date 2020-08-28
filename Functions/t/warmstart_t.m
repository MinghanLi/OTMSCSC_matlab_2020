function model = warmstart_t(X_Fold,model)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
disp('Warm-start...');
 %% prealignment
disp('preAlignment...');
[height,wide,n] = size(X_Fold);
TX = X_Fold; Tau = zeros(model.type,n);
[TX,Tau] = preAlign(TX,median(TX,3),Tau);

disp('initializing U...');
X = reshape(TX,[height*wide,n]);
%% PCA and L1MF for warmstart
[~, U, V] = inexact_alm_rpca(X);
r = model.r;
U = gather(U(:,1:r));
V = gather(V(:,1:r));
L = U*V';


%% Calculating the model
noi=X-L;
model.Sigma =var(noi(:))*50*0.05; %a simple initialization of mog parameters
model.D = InitDictNN(reshape(noi, [height,wide,n]), model.f_size);     % [f_size,f_size,K]

%% Calculating A and B
A=zeros(r,r,height*wide);B=zeros(r,height*wide);
for i=1:height*wide
    A(:,:,i)=(V'*V+0.001*eye(r))^-1*0.02;
    B(:,i)=V'* (L(i,:))'/0.02;
end
model.A=A; model.B=B;
model.U=U; model.v=V(1:2,:);model.r=r;
model.FA=zeros(length(model.f_size));
model.FB=zeros(max(model.f_size)^2,length(model.f_size));
disp('Warm start is over. ');
end

function  Filters  = InitDictNN( R, f_size )
% pca ≥ı º—°‘ÒD
FInd = FilterInd(f_size);
maxf_size = max(f_size);
K = length(f_size);
Patches = Video2Patch( R, maxf_size );
Patches = Patches- mean(Patches,2); 
temp = Patches*Patches';
[U, ~, ~] = svd(temp);
D = Normalize(max(U(:,2:K+1),0)); 
Filters = zeros(size(D)); Filters(FInd) = D(FInd);
end