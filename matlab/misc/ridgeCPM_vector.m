function [q_s, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total] = ridgeCPM_vector(all_mats, all_vectors, all_behav, thresh, famid, v_alpha, lambda, k, seed)
    %ridgeCPM Connectome-based predictive modeling using univariate
    %feature selection and ridge regression, additional handling
    %vector-shaped data
    %
    %   [q_s, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total] = ridgeCPM_vector(all_mats, all_vectors, all_behav, thresh, v_alpha, lambda, k, seed)
    %
    %   Input:      all_mats,           connectome of all the subjects and tasks
    %                                   [regions x regions x subjects x tasks]
    %
    %               all_vectors,         additional features for each subject
    %                                   [n_features x subjects]
    %
    %               all_behav,          behavior of all the subjects
    %                                   [subjects x 1]
    %
    %               thresh,             p-value threshold used for univariate
    %                                   feature selection
    %
    %               v_alpha(optional),  value of the alpha parameter in elastic 
    %                                   net, default is 1e-6 which makes the
    %                                   regression method to be ridge
    %                                   regression, v_alpha=1 makes it lasso.
    %
    %               lambda(optional),   value of the lambda, if not provided, 
    %                                   cross-validation will be used
    %
    %               k(optional),        number of folds in k-fold cross
    %                                   validation, default is 10
    %
    %               seed(optional),     random seed, default is 665
    %
    %   Output:     q_s,                cross-validated R^2 between predicted
    %                                   value with ground truth (all_behav)
    %
    %               r_pearson,          WRONG! direct pearson correlation 
    %                                   between predicted value with ground
    %                                   truth, only kept here for comparison
    %                                   will be removed afterwards
    %
    %               r_rank,             cross-validated spearman correlation
    %                                   between predicted value with ground 
    %                                   truth (all_behav)
    %
    %               y,                  predicted value for all the subjects
    % 
    %               coef_total,         regression coefficients of all the edges
    %                                   in all the k folds
    %
    %               coef0_total,        regression intercept in all the k folds
    %
    %               lambda_total,       penalty parameter chosen at each
    %                                   iteration
    %       
    %   Siyuan Gao, Yale University, 2018-2019
    %   Adapted by Abby Greene, 2019
    
    %% initialization
    if nargin < 3
        error('not enough arguments, please check the help')
    end
    
    if ~exist('k', 'var')
        k = 10;
    end
    
    if ~exist('seed', 'var')
        seed = 665;
    end
    
    if ~exist('v_alpha', 'var') || length(v_alpha) ~= 1 
        v_alpha = 1e-9;
    end
    
    % set the number of lambda to search
    numLambda = 16;
    
    % parse dimensions of the data
    if ~isempty(all_mats)
        num_sub_total = size(all_mats, 3);
        num_node = size(all_mats, 1);
        num_task = size(all_mats, 4);
    else % if only vector input, set variables based on size of vector
        num_sub_total = size(all_vectors,2);
        num_node = size(all_vectors,1);
        num_task = size(all_vectors,3);
    end
    num_feature = size(all_vectors, 1);
    
    % determine number of edges based on whether matrix is symmetric
    is_sym = issymmetric(all_mats(:, :, 1, 1));
    if is_sym && ~isempty(all_mats)
        num_edge = num_node * (num_node - 1) / 2;
    elseif ~is_sym && ~isempty(all_mats)
        num_edge = num_node * num_node;
    else % if only vector input
        num_edge = 0;
    end

    coef_total = zeros(num_edge*num_task+num_feature, k); %store all the coefficients
    coef0_total = zeros(1, k); % store all the intercept
    lambda_total = zeros(1, k); % store all the lambda

    %% convert connectivity to edge matrix (could made easier by squareform) - only if there's an all_mats input; if only vector, skip this
    if ~isempty(all_mats)
        all_edges = zeros(num_edge, num_sub_total, num_task);
        for i_sub = 1 : num_sub_total
            for j_task = 1 : num_task
                if is_sym
                    all_edges(:, i_sub, j_task) = squareform(tril(all_mats(:, :, i_sub, j_task), -1));
                else
                    all_edges(:, i_sub, j_task) = reshape(all_mats(:, :, i_sub, j_task), [], 1);
                end
            end
        end
        all_edges = permute(all_edges, [1, 3, 2]);
        all_edges = reshape(all_edges, [], num_sub_total);
    else
    end
    
    %% main
    
    y = zeros(num_sub_total, 1);
    rng(seed);
%     indices = crossvalind('Kfold', num_sub_total, k);
    indices = kfold_family(famid,k,seed); % ASG added this to run kfold_family idx generation on the fly
    
    for i_fold = 1 : k
        tic
        fprintf('%dth fold\n', i_fold);
        
        test_idx = (indices==i_fold);
        train_idx = (indices~=i_fold);
        
        if ~exist('all_edges','var') % if only vector input
            train_mats = all_vectors(:, train_idx);
            test_mats = all_vectors(:, test_idx);
        elseif isempty(all_vectors) % if only matrix input
            train_mats = all_edges(:, train_idx);
            test_mats = all_edges(:, test_idx);
        else % if both vector and matrix input
            train_mats = [all_edges(:, train_idx); all_vectors(:, train_idx)];
            test_mats = [all_edges(:, test_idx); all_vectors(:, test_idx)];
        end
        
        train_behav = all_behav;
        train_behav(test_idx) = [];
        
        % first step univariate edge selection
        [~, edge_p] = corr(train_mats', train_behav);
        edge_idx = find(edge_p < thresh);
        disp(['edge size = ' num2str(size(edge_idx))]);
     
        % build model on TRAIN subs
        % Choosing lambda
        if ~exist('lambda', 'var') || numel(lambda) < 1
            % automatically choose lambda
            lmax = computeLambdaMax(train_mats(edge_idx, :)', train_behav, [],...
                0.01, true);
            lmin = lmax * 0.001;
            loghi = log(lmax);
            loglo = log(lmin);
            lambda_squence = exp(linspace(loghi,loglo,numLambda));
            opt_lambda = find_lambda(train_mats, train_behav,...
                lambda_squence, thresh, 5); % TIME-CONSUMING
        elseif numel(lambda) > 1 % user-specified lambda sequence
            opt_lambda = find_lambda(train_mats, train_behav,...
                lambda, thresh);
        else % user-specified lambda value
            opt_lambda = lambda;
        end
        lambda_total(i_fold) = opt_lambda; % record the chosen lambda
        
        % Second step train the ridge model with the optimal lambda
        [beta, fit_info] = lasso(train_mats(edge_idx, :)',...
            train_behav, 'Alpha',v_alpha, 'Lambda', opt_lambda);
        intercept = fit_info.Intercept;

        % run model on TEST sub with the best lambda parameter
        
        y(test_idx) = test_mats(edge_idx, :)'*beta+intercept;
        
        coef_total(edge_idx, i_fold) = beta;
        coef0_total(:, i_fold) = intercept;
        toc
    end
    
    % compare predicted and observed behaviors
    [r_pearson, ~] = corr(y, all_behav);
    [r_rank, ~] = corr(y, all_behav, 'type', 'spearman');
    mse = sum((y - all_behav).^2) / num_sub_total;
    q_s = 1 - mse / var(all_behav, 1);
end

