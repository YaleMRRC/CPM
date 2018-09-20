function [r,p,pmask,mdl]=cpm_train(x,y,pthresh)
% call from predict_behavior

% feature selection from all train subs
[r,p]=corr(x',y);
pmask=(+(r>0))-(+(r<0)); % BINARIZED
pmask=pmask.*(+(p<pthresh));

% get summary features for each train sub
for i=1:size(x,2)
    % TODO: options for sum v mean, pos v neg v all
    summary_feature(i)=nanmean(x(pmask>0,i))-nanmean(x(pmask<0,i)); % 1. extract features
%     summary_feature(i)=nanmean(x(:,i).*+(pmask>0))-nanmean(x(:,i).*+(pmask<0)); % 2. keep zeros
end

% fit behavior to summary features
% mdl=polyfit(summary_feature,y',1); % 1. polyval
mdl=robustfit(summary_feature,y'); % 2. robust regression
    
    