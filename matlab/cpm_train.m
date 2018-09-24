function [r,p,pmask,mdl]=cpm_train(x,y,pthresh)
% call from cpm_cv
% Train a model using CPM framework
% Accepts the explanatory variable, "x", the target variable "y" and the
% significance threshold "pthresh"

% Feature selection from all train subs
[r,p]=corr(x',y);
% Binarized mask of features
pmask=(+(r>0))-(+(r<0)); 
pmask=pmask.*(+(p<pthresh));

% Get summary features for each train sub
for i=1:size(x,2)
    summary_feature(i)=nanmean(x(pmask>0,i))-nanmean(x(pmask<0,i));
end

% Fit behavior to summary features
mdl=robustfit(summary_feature,y');    
    