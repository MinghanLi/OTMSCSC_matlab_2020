function [Filters, model]=update_filter(x,M_ori,model,FInd)
%% updata auxiliary variable A and B for U
A=model.FA;  % [3p^2, 3p^2]
B=model.FB;  % [3p^2,1]
frame = model.frame;
ro = (frame-1)/frame;
h=1/frame;
f_size=model.f_size; f_maxsize=max(f_size);
M_pad = padarray(M_ori,[1 1 0]*((f_maxsize)-1)/2,'symmetric','both');
% extract f_maxsize x f_maxsize patches
tic
for i=1:length(f_size)
    M(:,(i-1)*f_maxsize^2+1:i*f_maxsize^2)=im2col(M_pad(:,:,i),[f_maxsize f_maxsize],'sliding')'; % [n, 3*p^2]
end
% M=im2col(M_pad,[f_maxsize f_maxsize],'sliding')'; % [n, 3*p^2]
toc
new_A=ro*A+h*(M'*M);
new_B=ro*B+h*M'*x;
D= reshape(new_A\new_B,[f_maxsize^2,length(f_size)]);
D = Normalize(max(D,0));

%% tensor coding for upspeeding
% temp= 1/ro*A*M'; % [3*p^2, n]
% new_A=1/ro*A-(ro+h*temp*M)\(h*temp*M*A);
% new_B=B*ro+h*M'*x;
% D=reshape(new_A*new_B,[f_maxsize^2,length(f_size)]);
Temp = zeros(size(D)); Temp(FInd) = D(FInd);
Filters = reshape(Temp,[f_maxsize,f_maxsize,length(f_size)]);
imshow([Filters(:,:,1)*10,Filters(:,:,2)*10,Filters(:,:,3)*10])
%% output 
model.FA=new_A;
model.FB=new_B;