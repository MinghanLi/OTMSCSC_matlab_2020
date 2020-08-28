function Ind = FilterInd(f_size)
K = length(f_size); maxf_size = max(f_size);
ind = (maxf_size-f_size)/2;
IndMax = zeros(maxf_size,maxf_size,K);
for k=1:K
    IndMax(1+ind(k):end-ind(k),1+ind(k):end-ind(k),k) = 1;
end
Ind = IndMax==1; 
end