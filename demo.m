clc; clear; close all
currentFolder = pwd; addpath(genpath(currentFolder))

%% initial mask param
param.alpha = 1; % alpha controls the thresholding, smaller lambda means more mask
param.rho = 1.5;
param.istransform=0;
param.weight=1;
model.use_gpu = 'Ture';
model.f_size =[13 9 5];     
model.r=1; model.MaxIter=1; model.isMask=1;
model.preAlignment = 0; model.type=3; % 3:rigid, 4:simi, 6:affine, 8:projective
model.b=0.001*ones(1,length(model.f_size)); 
model.Flamda=0.01;
model.display =1 ;
model.interval_alignB=2000; % 0 means no align B, larger have faster speed, but may degrade perfomance
model.strNum=50;
model.derain=2; % 1 means Derain=H*B+(1-H)*F, 2 means Derain=X-R

model.video_path='/home/lmh/Downloads/OTMSCSC/OTMSCSC_matlab/input/CDNet/bad_weather/skating/';

%add your video image path 
video_path = [model.video_path,'input/'];
[t, sigma, b] = run_video(video_path,param,model);
plot(t)
% for i=1:300

% t1(i)=sum(t(1:i));
% end
% xlswrite('F:\NoiseModel\Background\Derain\DerainTPAMI\time.xls',t,'t_london')
% load 4080GT
% metric = metrics(double(groundtruth),Derain*255)
