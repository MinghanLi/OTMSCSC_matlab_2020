clc; clear; close all
currentFolder = pwd; addpath(genpath(currentFolder))

%% initial mask param
% alpha controls the thresholding, smaller lambda means more mask
param.alpha = 1; 
param.rho = 1.5;
param.weight=1;
% Although we also support CPU, we recommend you to run the code in devices with GPU. 
model.use_gpu = 'Ture';
model.batch_size = 2;
model.f_size =[13 9 5];  % filter size
model.r=1; model.MaxIter=1; 
model.preAlignment = 0; 

% if there are not moving objects in videos, you can assign isMask=0 to obtain faster speed.
% Of course, isMask=1 can deal with all videos.
model.isMask=1;
% please assign 1 for transformed videos
param.istransform = 0;
% 3:rigid, 4:simi, 6:affine, 8:projective
model.type=3; 
model.b=0.001*ones(1,length(model.f_size)); 
model.Flamda=0.01;
% 0 means no align B, larger have faster speed, but may degrade perfomance
model.interval_alignB=2000; 
model.strNum=50;
model.derain=2; % 1 means Derain=H*B+(1-H)*F, 2 means Derain=X-R
model.display=1;

%add your video image path
model.video_path='/home/lmh/Downloads/OTMSCSC/OTMSCSC_matlab_2020/input/skating/';
video_path = [model.video_path,'input/'];
[t, sigma, b] = run_video(video_path,param,model);

% load GroundTruth to calculate metrics (PSNR, SSIM)
% metric = metrics(double(groundtruth),Derain*255)
