function [y_predict]=cpm_test(x,mdl,pmask)
% do prediction using model built in cpm_train

% do prediction
for i=1:size(x,2)
    summary_feature(i)=nanmean(x(pmask>0,i))-nanmean(x(pmask<0,i)); % 1. extract features
%     summary_feature(i)=nanmean(x(:,i).*+(pmask>0))-nanmean(x(:,i).*+(pmask<0)); % 2. keep zeros
%     y_predict(i)=polyval(mdl,summary_feature(i)); % 1. polyval
    y_predict(i)=mdl(2)*summary_feature(i) + mdl(1); % 2. robust regression
end
