function  Patches  =  Im2Patch( D,PatSize )
TotalPatNum = (size(D,1)-PatSize+1)*(size(D,2)-PatSize+1);                  %Total Patch Number in the image
Patches   =   zeros(PatSize*PatSize, TotalPatNum,'single');                      %Current Patches
%Patches   =  zeros(PatSize*PatSize, TotalPatNum);   
k   =   0;
for j  = 1:PatSize
    for i  = 1:PatSize
        k           =  k+1;
        temp        =  D(i:end-PatSize+i,j:end-PatSize+j);  % 按列展开[PatSize*PatSize, TotalPatNum]     
        Patches(k,:)=  temp(:)';
    end
end
 