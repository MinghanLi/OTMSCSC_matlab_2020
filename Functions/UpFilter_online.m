function [Filters, model]=UpFilter_online(X_ori, model, FInd)
X_ori = X_ori-min(X_ori(:)); f_size=model.f_size; f_maxsize=max(f_size);
X_pad = padarray(X_ori,[1 1 0]*((f_maxsize)-1)/2,'symmetric','both');
% extract f_maxsize x f_maxsize patches
X=im2col(X_pad,[f_maxsize f_maxsize],'sliding');
X=X-repmat(mean(X),[size(X,1) 1]);
X=X ./ sqrt(sum(X.^2));

param.K=length(f_size);  % learns a dictionary with 100 elements
param.lambda=0.15;
param.numThreads=4; % number of threads
param.batchsize=256;
param.D=model.D;
param.approx=0; 
param.iter=50;  % let us see what happens after 1000 iterations.

%%%%%%%%%% FIRST EXPERIMENT %%%%%%%%%%%
Fmodel.A=model.FA;
Fmodel.B=model.FB;
Fmodel.iter=100;
[D,Fmodel]= mexTrainDL(X,param,Fmodel);
D = Normalize(max(D,0));
model.D = D; model.FA=Fmodel.A; model.FB=Fmodel.B;
Temp = zeros(size(D)); Temp(FInd) = D(FInd);
Filters = reshape(Temp,[f_maxsize,f_maxsize,length(f_size)]);


