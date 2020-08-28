
video_path='/home/lmh/Downloads/OTMSCSC/data_result/spac/synthetic/';
name_gt='GT/a3_GT/';
name_derain = '3Ren/a3_ren/'; 

path_gt=dir([video_path, name_gt]);
path_derain = dir([video_path, name_derain]);
len=length(path_derain)-2;
first_frame =  im2double(imread([video_path,name_gt,path_gt(3).name]));
img_size = size(first_frame);

GT = zeros(img_size(1), img_size(2),img_size(3),  len);
Derain = zeros(img_size(1), img_size(2), img_size(3), len);
metric = zeros(len, 6);
for i = 1:len
    GT= im2double(imread([video_path,name_gt,path_gt(i+2).name]));
    Derain =  im2double(imread([video_path,name_derain,path_derain(i+2).name]));
    metric(i, :) = metrics(GT,Derain);
end
sum(metric, 1)/len
    

