function [mdl_summary, mdl_features, y_test, y_predict]=cpm_cv(x,y,pthresh,kfolds)
% call from predict_behavior
% can provide previously defined models (mdl1,mdl2) as input
% can convert binarized to weighted--see 3 chunks below

% TODO: catch kfold outside
% TODO: save features

% splitting data
nsubs=size(x,2);
nfeats=size(x,1);
randinds=randperm(nsubs);
ksample=floor(nsubs/kfolds);

% maybe preallocate

fprintf('\n# Running over %1.0f Folds. Performing fold # ',kfolds);
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

mdl_summary=0; %TBD
mdl_features=0;%TBD

fprintf('\n');
    
    
    
