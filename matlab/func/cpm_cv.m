function [y_predict]=cpm_cv(x,y,pthresh,kfolds)
% Runs cross validation for CPM
% x            Predictor variable
% y            Outcome variable
% pthresh      p-value threshold for feature selection
% kfolds       Number of partitions for dividing the sample
% y_test       y data used for testing
% y_predict    Predictions of y data used for testing

% Split data
nsubs=size(x,2);
randinds=randperm(nsubs);
ksample=floor(nsubs/kfolds);

y_predict = zeros(nsubs, 1);
% Run CPM over all folds
fprintf('\n# Running over %1.0f Folds.\nPerforming fold no. ',kfolds);
for leftout = 1:kfolds
    fprintf('%1.0f ',leftout);
    
    if kfolds == nsubs % doing leave-one-out
        testinds=randinds(leftout);
        traininds=setdiff(randinds,testinds);
    else
        si=1+((leftout-1)*ksample);
        fi=si+ksample-1;
        
        testinds=randinds(si:fi);
        traininds=setdiff(randinds,testinds);
    end
    
    % Assign x and y data to train and test groups 
    x_train = x(:,traininds);
    y_train = y(traininds);
    x_test = x(:,testinds);
    
    % Train Connectome-based Predictive Model
    [~, ~, pmask, mdl] = cpm_train(x_train, y_train,pthresh);
    
    % Test Connectome-based Predictive Model
    [y_predict(testinds)] = cpm_test(x_test,mdl,pmask);
end
