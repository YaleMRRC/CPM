function [x, y] = cpm_check_errors(x, y, folds)
% Checks that input data is in a format usable by CPM

% Check that x data are in the required format - 2D or 3D
if (ndims(x)~=2) && (ndims(x)~=3)
   error('Data should have two or three dimensions') 
end

% Check that x data contain more than one element
if size(x,1)==1
    error('Single feature detected.')
end

% Check that there the same number of subjects in x as in y
if size(y,1)==1 % If behavioral scores are row vector, reformat to be column vector
    y = y';
end

if size(x,ndims(x))~=size(y,1)
    error('There are NOT the same number of subjects in the data and behavior variable')
end

% Check to make sure there are at least ten subjects in the input data
if size(x,ndims(x))<10
    warning('The CPM code requires >10 subjects to function properly; sound results likely require >>10.')
end

% Check to make sure you have more subjects than folds
if size(x,ndims(x))<folds
    warning('You must have more subjects than folds in your cross validation. Please check the help documentation.')
end

% Check whether x is symmetric across first two dimensions
if ndims(x)==3
    if size(x,1)~=size(x,2)
        warning('Please make sure, if intended, that data is an NxN connectivity matrix')
    end
end

% Check for nodes with values of 0 (missing nodes within a subject)
row_sum = squeeze(sum(abs(x), 2));
zero_node = sum(row_sum==0);
zero_node_subjects = sum( zero_node>0);
if zero_node_subjects>0
    warning('Data: %d subjects have missing nodes. Please check your data.',zero_node_subjects)
end

% Check for Inf or NaN 
if length(find(isinf(x)))>0
    warning('You have Inf values in your matrices. Please check your data.')
end

if length(find(isnan(x)))>0
    warning('You have NaNs in your matrices. Please check your data.')
end

% If data are 3D, convert to 2D
% If matrix is symmetric, only upper triangle is taken
if ndims(x)==3
    if issymmetric(x(:,:,1))
        s=size(x,1);
        for i=1:size(x,3)
            data_tmp=x(:,:,i);
            m(:,i)=data_tmp(logical(triu(ones(s,s),1)));
        end
    else
        m=reshape(x,size(x,1)*size(x,2),size(x,3));
    end
    x=m;
end

end

