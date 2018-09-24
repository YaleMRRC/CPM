function [y_predict_reshape]=cpm_cv(x,y,pthresh,kfolds)
% Runs cross validation for CPM
% x            Predictor variable
% y            Outcome variable
% pthresh      p-value threshold for feature selection
% kfolds       Number of partitions for dividing the sample
% y_test       y data used for testing
% y_predict    Predictions of y data used for testing

% Split data
nsubs=size(x,2);
nfeats=size(x,1);
randinds=randperm(nsubs);
ksample=floor(nsubs/kfolds);

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
    nsubs_in_fold=length(testinds);
    
    % Assign x and y data to train and test groups 
    x_train = x(:,traininds);
    y_train = y(traininds);
    x_test = x(:,testinds);
    y_test(leftout,1:nsubs_in_fold) = y(testinds);
    
    % Train Connectome-based Predictive Model
    [r,p,pmask,mdl] = cpm_train(x_train, y_train,pthresh);
    
    % Test Connectome-based Predictive Model
    [y_predict(leftout,1:nsubs_in_fold)]=cpm_test(x_test,mdl,pmask);
end

y_test_reshape(randinds)=reshape(y_test',[],1);
y_predict_reshape(randinds)=reshape(y_predict',[],1);
