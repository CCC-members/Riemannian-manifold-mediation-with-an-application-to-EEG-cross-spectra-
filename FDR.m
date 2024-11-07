function [pID,pN] = FDR(p,q)
% FORMAT [pID,pN] = FDR(p,q)
% 
% p   - vector of p-values
% q   - False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold
%______________________________________________________________________________
% $Id: FDR.m,v 2.1 2010/08/05 14:34:19 nichols Exp $


p = p(isfinite(p));  % Toss NaN's
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*q/cVID)));
if isempty(pID), pID=0; end
pN = p(max(find(p<=I/V*q/cVN)));
if isempty(pN), pN=0; end

% 
% The following function will supply p-value thresholds which control the expected FDR at a specified rate.
% 
% FDR.m 
% This function takes a vector of p-values and a FDR rate. It returns two p-value thresholds, one based on an assumption of independence or positive dependence, and one that makes no assumptions about how the tests are correlated. For imaging data, an assumption of positive dependence is reasonable, so it should be OK to use the first (more sensitive) threshold.
% Here is one example of how to use this function with a T image.
% 
%  
%       Timg     = % read t image, with df degrees of freedom
%       Timg(~isfinite(Timg(:))) = [];
%  
%       P        = 1-cdf('T',Timg(:),df);
%  
%       [pID pN] = FDR(P,0.05);
%  
%       tID      = icdf('T',1-pID,df)     % T threshold, indep or pos. correl.
%       tN       = icdf('T',1-pN,df)      % T threshold, no correl. assumptions