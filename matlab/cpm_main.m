function [mdl_summary, mdl_features, y_predict, performance] = cpm_main(x,y,varargin)
% Performs CPM
% x is n x nsubs, or n x m x nsubs
% y is i x nsubs
% returns m, a structure containing the predictive edges and associated p and r vals
% e.g., [m1,m2,ypred,perf]=cpm_main(data,gF,'pthresh',0.05,'kfolds',2);

%% Parse input
p=inputParser;
defaultpthresh=0.05;
defaultkfolds=2; % TODO: loop over vec

addRequired(p,'data',@isnumeric); % must be n x n x subs or nfeatures x nsubs
addRequired(p,'behavior_data',@isnumeric); % must be n x nsubs
addParameter(p,'pthresh',defaultpthresh,@isnumeric);
addParameter(p,'kfolds',defaultkfolds,@isnumeric);

parse(p,x,y,varargin{:});

clearvars -except x y p

pthresh = p.Results.pthresh;
kfolds = p.Results.kfolds;


%% Errors

if size(x,ndims(x)) ~= size(y,1)
    warning('warning TBD')
end

% TODO: stop if ndims(x) > 3d
% TODO: catch if only one feature - WHYYYY

%% Make data 2D

if ndims(x)~=2 % 3d
    if issymmetric(x(:,:,1))
        s=size(x,1);
        for i=1:size(x,3)
            data_tmp=x(:,:,i);
            masked_data(:,i)=data_tmp(logical(tril(ones(s,s),-1)));
        end
    else
        masked_data=reshape(x,size(x,1)*size(x,2),size(x,3));

    end
end


%% Run train / test
[mdl_summary,mdl_features,y_test, y_predict]=cpm_cv(masked_data,y,pthresh,kfolds);
% TODO: choose + return mdl_summary, mdl_features

%% Check performance
performance=corr(y_predict(:),y_test(:));

end
