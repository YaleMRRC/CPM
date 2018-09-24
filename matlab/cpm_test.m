function [y_predict]=cpm_test(x,mdl,pmask)
% Test a Connectome-based Predictive Model using previously trained model
% x            Predictor variable
% mdl          Coefficient fits for linear model relating summary features to y
% pmask        Mask for significant features
% y_predict    Predicted y values

% For each subject, create summary feature and use model to predict y
for i=1:size(x,2)
    summary_feature(i)=nanmean(x(pmask>0,i))-nanmean(x(pmask<0,i));
    y_predict(i)=mdl(2)*summary_feature(i) + mdl(1); 
end
