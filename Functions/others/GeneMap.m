function GMap = GeneMap(Map)
[h,w,K] = size(Map);
for k=1:K
    OutZK = reshape(Map(:,:,k),h*w,1);
    ind = OutZK<=00.005; 
    ZeroRate = sum(ind)/length(OutZK); 
    % subplot(1,2,1)
    % hist(OutZK(~ind),100)
    mu = 0; % mu = mean(OutZK(index));                                       % find the distribution expected value
    sigma = sum(abs(OutZK(~ind)-mu))/sum(~ind);                            % find the distribution std
    
    % Generate Laplacian noise
    ind = rand(h*w,1) > min(ZeroRate+0.15*(K+1-k),0.99);                                                % find nonzero location
    u = rand(sum(ind),1) -0.5;
    yrgb = zeros(h*w,1);
    temp = abs(mu - sigma * sign(u).* log(1- 2* abs(u)));
    % subplot(1,2,2)
    % hist(temp,100)
    yrgb(ind) = temp; yrgb(find(ind==1)+1) = temp;
    % yrgb(find(ind==1)+2) = temp; yrgb(find(ind==1)+h) = temp;
    GMap(:,:,k) = reshape(yrgb(1:h*w),[h,w]);
end
% for i=1:rgb % rgb channels
%     for k = 1:K
%         GeneZpat = Video2Patch( reshape(TempGeneOutZ(:,:,i,:,k),[h+f_size-1,w+f_size-1,n]),f_size ); 
%         GeneOutRS(:,:,i,:,k) = reshape(Filters((k-1)*f_size^2+1:k*f_size^2,i)'* GeneZpat,[h,w,n]);
%     end
% end
% GeneOutRS = max(GeneOutRS-0.001,0);
end
