function [x, y] = cpm_check_errors(x, y, folds)
%This function checks the dataset to make sure it is in a format usable by
%CPM


%INPUTS
%           data                        A 2D or 3D matrix, N x nSubjects
%                                       OR N x M x nSubjects
%
%           behavior                    2D column vector of behavioral scores,
%                                       Subjects x 1 (each row is a
%                                       subject's behavioral score)
%

%OUTPUT
%
%           pass                        This is a 1 (pass) or 0 (fail)
%           Message                     This prints a message


%1) is data in the required format? 2D or 3D x, column vector

if (ndims(x)~=2) && (ndims(x)~=3)
   error('Data should have two or three dimensions') 
end

if (size(x,1)==1 && size(x,1)==1)
    error('Single feature detected.')
end

%2) Number of subjects in 3D arrays is same as the behavioral score
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

%3) Check there are the same number of nodes, that x is symmetric
%across first TWO dimensions

if size(x,1)~=size(x,2)
    warning('Please make sure, if intended, that data is an NxN connectivity matrix')
end

%4) Check for nodes with values of 0 (missing nodes within a subject)
row_sum = squeeze(sum(abs(x), 2));
zero_node = sum(row_sum==0);
zero_node_subjects = sum( zero_node>0);
if zero_node_subjects>0
    warning('Data: %d subjects have missing nodes. Please check your data.',zero_node_subjects)
end

%4) Check for Inf or NaN 
if length(find(isinf(x)))>0
    warning('You have inf values in your matrices. Please check your data.')
end

if length(find(isnan(x)))>0
    warning('You have NaNs in your matrices. Please check your data.')
end

%5) Zero out the diagonal; maybe do this later?

% Maybe check CV parameters to make sure we have enough subjects


%% Make data 2D

if ndims(x)~=2 % 3d
    if issymmetric(x(:,:,1))
        s=size(x,1);
        for i=1:size(x,3)
            data_tmp=x(:,:,i);
%             m(:,i)=data_tmp(logical(tril(ones(s,s),-1)));
            m(:,i)=data_tmp(logical(triu(ones(s,s),1)));
        end
    else
        m=reshape(x,size(x,1)*size(x,2),size(x,3));
    end
    x=m;
end

end


