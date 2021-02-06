clc,clear all
% type = 1 means to calculate all videos in a folder
% type =0 means to only calcuate a specific video
type=0;

if type==1
    videos_derain_path='/home/lmh/Downloads/OTMSCSC/dataset/YTVOS2019_all_frames/9SLDNet/';
    videos_gt_path='/home/lmh/Downloads/OTMSCSC/dataset/YTVOS2019_all_frames/GT/';
    videos_dir = dir(videos_derain_path);
    videos_len = length(videos_dir)-2;

    metric = zeros(videos_len, 2);
    for i = 1:videos_len
        cur_gt_path = dir([videos_gt_path, videos_dir(i+2).name]);
        cur_derain_path = dir([videos_derain_path, videos_dir(i+2).name]);
        cur_derain_len = length(cur_derain_path) - 2;
        cur_metric = zeros(cur_derain_len, 2);
        temp = i
        a=cur_gt_path(1+2+4).name
        b=cur_derain_path(1+2).name
        for j = 1:cur_derain_len
            GT = im2double(imread([videos_gt_path, videos_dir(i+2).name, '/', cur_gt_path(j+2+4).name]));
            Derain =  im2double(imread([videos_derain_path,videos_dir(i+2).name, '/',cur_derain_path(j+2).name]));
            [h, w, ch] = size(Derain);
            GT = imresize(GT, [h, w]);
            cur_metric(j, :) = metrics(GT,Derain);
        end
        metric(i, :) = sum(cur_metric, 1) / cur_derain_len;
    end
    sum(metric, 1)/videos_len
else
    video_derain_path='/home/lmh/Downloads/OTMSCSC/data_result/spac/synthetic/1Garg/';
    video_gt_path='/home/lmh/Downloads/OTMSCSC/data_result/spac/synthetic/GT/';
    folder='a4';

    path_gt=dir([video_gt_path, folder, '/']);
    path_derain = dir([video_derain_path, folder, '/']);
    len=length(path_derain)-4;
    % len = 60;

    metric = zeros(len, 2);
    output_save = zeros(len, 3);
    for i = 1:len
        GT= im2double(imread([video_gt_path,folder, '/', path_gt(i+2).name]));
        Derain =  im2double(imread([video_derain_path,folder, '/', path_derain(i+2).name]));
        % [h, w, ch] = size(Derain);
        % GT = imresize(GT, [h, w]);
        metric(i, :) = metrics(GT,Derain);
    end
    sum(metric, 1)/len
    output_save(:, 2:3) = metric;
    for i = 1:len-1
        output_save(i,1) = i;
    end
    save_path = [video_derain_path, folder, '.mat'];
    save(save_path, 'output_save')
end

