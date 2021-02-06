% load Tau
% load r_ali
% for i=1:400
%     I=rgb2ycbcr(input(:,:,:,i));
%     I(:,:,1)=warpImg(I(:,:,1),Tau(:,i));
%     t_input(:,:,:,i)=ycbcr2rgb(I);
% end
video_path='F:\NoiseModel\Background\Derain\MS-CSConline\Online-MS-CSC-Rain-Streak-Removal\input\jitter\t_tree\result\';
name='GT\';
% for i=1:size(rain,4)
%     rg(:,:,:,i)=repmat(rgb2gray(rain(:,:,:,i)),[1,1,1]);
% end
input=groundtruth;
% name='result\';
% input=OutDeRain;
for i =1:size(input,4)
    I = input(:,:,:,i);
    if i<10
        imwrite(I,[video_path,name,'000',int2str(i),'.bmp']);
    elseif i>=10 && i <100
        imwrite(I,[video_path,name,'00',int2str(i),'.bmp']);
    elseif i>=100 && i <1000
        imwrite(I,[video_path,name,'0',int2str(i),'.bmp']);
    else
         imwrite(I,[video_path,name,int2str(i),'.bmp']);
    end
end