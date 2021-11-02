function  [t, sigma, b] = run_video(video_path, param, model)
%% init
channel=1;  % transfer images from RGB to YCbCr color space, where Y=1,Cb=2,Cr=3;
path=dir(video_path);
len=length(path);
display=model.display;

%% warmstart 
firIm=imread([video_path,path(5).name]);
imgsize=size(firIm);


%% init mask
 param.Mask = true([imgsize(1),imgsize(2),2]);
 % param.MaskOut = false([imgsize(1),imgsize(2),2]);
 
%% main 
% transfer images from RGB to YCbCr color space
firIm = rgb2ycbcr(im2double(firIm));
model.lastframe = firIm(:,:,channel); 
for i=1:ceil((len-2)/model.batch_size)
    model.frame = i;
    
    if mod(i,10)==0||i==1
        disp(['Calculating the model of the ',num2str(i),'th frame']);
    end
    s=tic;
    %% warm up
    if param.istransform == 0 && i==1
        strNum=min(model.strNum, len-2);
        ind=sort(randperm(len-2, strNum)+i-1);
        X_warm=zeros(imgsize(1),imgsize(2),strNum);
        for j = 1:strNum
            im_warm=im2double(imread([video_path,path(ind(j)+2).name]));
            im_warm_Ycbcr = rgb2ycbcr(im_warm);
            X_warm(:,:,j)=im_warm_Ycbcr(:,:,channel);
        end
        model = warmstart(X_warm,model);
    elseif param.istransform == 1  && mod(i-1,model.interval_alignB)==0
        strNum=model.strNum;
        ind=sort(randi([max(i-strNum, 1), min(i+strNum, len)], strNum, 1));
        X_warm=zeros(imgsize(1),imgsize(2),strNum);
        for j = 1:strNum
            im_warm=im2double(imread([video_path,path(ind(j)+2).name]));
            im_warm_Ycbcr = rgb2ycbcr(im_warm);
            X_warm(:,:,j)=im_warm_Ycbcr(:,:,channel);
        end
        model = warmstart_trans(X_warm,model);
    end
    clear X_warm

    n_frames = min((i+1)*model.batch_size, len-2) - i*model.batch_size;
    im_cur_Ycbcr = zeros([imgsize, n_frames]);  % [h, w, 3, n]
    for j = 1:n_frames
        im_cur = im2double(imread([video_path,path(i*model.batch_size+2+j).name]));
        im_cur_Ycbcr(:, :, :, j) = rgb2ycbcr(im_cur); 
    end
    X = im_cur_Ycbcr(:, :, channel, :);
    
    % main function
    if param.istransform == 0
        [model,DerainY, Rain, Rains, omega] = OnlineMSCSC_new(X, param, model);
    else
        % align background via every interval frames
        [model,param,DerainY, Rain, Rains, omega] = OnlineMSCSC_trans(X, param, model);
    end
    im_cur_Ycbcr(:,:,channel,:) = DerainY;
    Derain = ycbcr2rgb(im_cur_Ycbcr);
    t(i)=toc(s);
    sigma(i)=model.Sigma;
    b(i,:)=model.b;
    
    %% display results
    if display == 1
        for j = 1:n_frames
            I = Derain(:,:,:,j);
            R = [Rain,Rains(:,:,1, j),Rains(:,:,2, j),Rains(:,:,3, j)]*30; imshow(R)
            % imshow([im_cur,Derain,repmat(omega,[1,1,3])]); 
            if i<10
                imwrite(I,[model.video_path,'result/000',int2str(i),'.bmp']);
                imwrite(R,[model.video_path,'rain_layer/000',int2str(i),'.bmp']);
            elseif i>=10 && i <100
                imwrite(I,[model.video_path,'result/00',int2str(i),'.bmp']);
                imwrite(R,[model.video_path,'rain_layer/00',int2str(i),'.bmp']);
            elseif i>=100 && i <1000
                imwrite(I,[model.video_path,'result/0',int2str(i),'.bmp']);
                imwrite(R,[model.video_path,'rain_layer/0',int2str(i),'.bmp']);
            else
                imwrite(I,[model.video_path,'result/',int2str(i),'.bmp']);
                imwrite(R,[model.video_path,'rain_layer/',int2str(i),'.bmp']);
            end
        end
     end
end



 
