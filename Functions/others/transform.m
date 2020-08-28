clc;clear 
load snow_tree 
load Tau
Tau=[Tau,Tau];
%input=data;
[h,w,~,n]=size(input);
for i=1:n
    Temp = rgb2ycbcr(input(:,:,:,i));
    Iwarp = warpImg(Temp(:,:,1),Tau(:,i));
    Temp(:,:,1)=uint8(Iwarp);
    t_input(:,:,:,i)=ycbcr2rgb(Temp);
end
implay(t_input)
save t_snow_tree.mat t_input