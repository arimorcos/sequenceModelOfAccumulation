function [out] = sortTimeMax(in)
%sortTimeMax.m function which sorts array n x t with n neruons and t
%timepoints according to time of peak activity 

[~,index] = max(in,[],2); %find time of max value
[~,index2] = sort(index,'descend'); %sort based on time of max activity 
out = in(index2,:);