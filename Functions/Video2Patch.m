function  Patches  =  Video2Patch( D, PatSize, use_gpu )
if(nargin<3) || isempty(use_gpu)
    use_gpu=1;
  end
TotalPatNum = (size(D,1)-PatSize+1)*(size(D,2)-PatSize+1)*size(D,3);                  %Total Patch Number in the image  
if use_gpu
    Patches   =   gpuArray(zeros(PatSize*PatSize, TotalPatNum,'single'));
    D = gpuArray(single(D));
else
    Patches   =  zeros(PatSize*PatSize, TotalPatNum); 
end
% video2patch
k   =   0;
for j  = 1:PatSize
    for i  = 1:PatSize
        k           =  k+1;
        temp        =  D(i:end-PatSize+i,j:end-PatSize+j,:);  % ����չ��[PatSize*PatSize, TotalPatNum]     
        Patches(k,:)=  temp(:)';
    end
end
% gpu2workspace
if use_gpu
    Patches = gather(Patches);
end
 