
n = 800;
output=zeros(360,400,3,n);
for i= 1:n
    move=floor(i/10);
    output(:,:,:,i)=input(:,1+move:400+move,:,i);
end
output=uint8(output);