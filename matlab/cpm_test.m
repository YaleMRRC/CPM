function [y_predict]=cpm_test(x,mdl,pmask)
% do prediction using model built in cpm_train

% do prediction
for i=1:size(x,2)
%     summary_feature(i)=mean(x(pmask,i))-mean(x(npmask,i));
      summary_feature(i)=mean(x(pmask>0,i))-mean(x(pmask<0,i));
    y_predict(i)=polyval(mdl,summary_feature(i));
end
