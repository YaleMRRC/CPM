function [y_test, y_predict]=cpm_cv(x,y,pthresh,kfolds)
% call from predict_behavior
% Runs cross validation based on CPM framework. Accepts explanatory
% variable "x", target variable "y", significance threshold "pthresh", and
% number of folds "kfolds" (2 =split half, 10 = ten fold etc.).


% splitting data
nsubs=size(x,2);
nfeats=size(x,1);
randinds=randperm(nsubs);
ksample=floor(nsubs/kfolds);

fprintf('\n# Running over %1.0f Folds.\nPerforming fold no. ',kfolds);

for leftout = 1:kfolds
    fprintf('%1.0f ',leftout);
    
    if kfolds == nsubs %loo
        testinds=randinds(leftout);
        traininds=setdiff(randinds,testinds);
    else % folds
        si=1+((leftout-1)*ksample);
        fi=si+ksample-1;
        
        testinds=randinds(si:fi);
        traininds=setdiff(randinds,testinds);
        % catch even split - use less - warning
    end
    nsubs_in_fold=length(testinds);
    
    % leave out subject from matrices and behavior
    
    x_train = x(:,traininds);
    y_train = y(traininds);
    x_test = x(:,testinds);
    y_test(leftout,1:nsubs_in_fold) = y(testinds);
    
    % train
    [r,p,pmask,mdl] = cpm_train(x_train, y_train,pthresh);
    
    %test
    [y_predict(leftout,1:nsubs_in_fold)]=cpm_test(x_test,mdl,pmask);
    
end

    
    
    
