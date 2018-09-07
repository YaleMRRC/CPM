%% Check that there are 3 input variables 
if nargin < 3
    error('There are not enough arguments; please check the help')
end

function [all_mats, all_behav] = check_data(all_mats, all_behav, folds)

%This function checks the dataset to make sure it is in a format usable by
%CPM


%INPUTS
%           all_mats                    A 3D matrix, Nodes x Nodes x Subjects
%
%           all_behav                   2D column vector of behavioral scores,
%                                       Subjects x 1 (each row is a
%                                       subject's behavioral score)
%

%OUTPUT
%
%           pass                        This is a 1 (pass) or 0 (fail)
%           Message                     This prints a message


%1) is there data in the required formats: 3D all_mats, column vector
%all_behavs

if length(size(all_mats))~=3
   warning('all_mats does NOT have three dimensions') 
end


%2) Number of subjects in 3D arrays is same as the behavioral score
if size(all_behav,1)==1 % If behavioral scores are row vector, reformat to be column vector
    all_behav = all_behav';
end

if size(all_mats,3)~=size(all_behav,1)
    warning('there are NOT the same number of subjects in the all_mats and all_behav variable')
end

% Check to make sure there are at least ten subjects in the input data
if size(all_mats,3)<10
    warning('the CPM code requires >10 subjects to function properly; sound results likely require >>10.')
end

% Check to make sure you have more subjects than folds
if size(all_mats,3)<folds
    warning('You must have more subjects than folds in your cross validation. Please check the help documentation.')
end

%3) Check there are the same number of nodes, that all_mats is symmetric
%across first TWO dimensions

if size(all_mats,1)~=size(all_mats,2)
    warning('please make sure all_mats is an NxN connectivity matrix')
end

%4) Check for nodes with values of 0 (missing nodes within a subject)
row_sum = squeeze(sum(abs(all_mats), 2));
zero_node = sum(row_sum==0);
zero_node_subjects = sum( zero_node>0);
if zero_node_subjects>0
    warning('all_mats: %d subjects have missing nodes. Please check your data.',zero_node_subjects)
end

%4) Check for Inf or NaN 
if length(find(isinf(all_mats)))>0
    warning('you have inf values in your matrices. Please check your data.')
end
if length(find(isnan(all_mats)))>0
    warning('You have NaNs in your matrices. Please check your data.')
end

%5) Zero out the diagonal; maybe do this later?

% Maybe check CV parameters to make sure we have enough subjects

end


