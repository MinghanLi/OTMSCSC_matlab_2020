function [ Mat ] = Normalize( Mat )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
dim = size(Mat,1);
Mat = Mat./repmat(sqrt(sum(Mat.^2))+eps,dim,1);
end

