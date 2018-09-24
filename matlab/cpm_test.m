function [y_predict]=cpm_test(x,mdl,pmask)
% call from cpm_cv
% Do prediction using model built in cpm_train
% Accepts pre-computed model, "mdl", the mask corresponding to the
% significant features, "pmask", and the explanatory variable "x"

% do prediction
for i=1:size(x,2)
    % 1. extract features
    summary_feature(i)=nanmean(x(pmask>0,i))-nanmean(x(pmask<0,i));
    % 2. robust regression
    y_predict(i)=mdl(2)*summary_feature(i) + mdl(1); 
end
