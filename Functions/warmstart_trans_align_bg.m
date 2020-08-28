function model = warmstart_trans_align_bg(X_Fold,model)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
disp('align background...');
 %% prealignment
% disp('preAlignment...');
[height,wide,n] = size(X_Fold);
TX = X_Fold; Tau = zeros(model.type,n);
[TX,~] = preAlign(TX,TX(:,:,ceil(n/2)),Tau);

% disp('initializing U...');
X = reshape(TX,[height*wide,n]);
%% PCA and L1MF for warmstart
[~, U, V] = inexact_alm_rpca(X);
r = model.r;
U = gather(U(:,1:r));
V = gather(V(:,1:r));
model.B_Fold = reshape(U*V(fix(n/2)+1,:)',[height, wide]);model.r=r;
disp('new align background is over. ');
end