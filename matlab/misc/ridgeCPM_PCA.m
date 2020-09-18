function [q_s, q_s_fold, r_pearson, r_rank, y, new_behav, pca_coef_total, ...
        coef_total, coef0_total, lambda_total] = ridgeCPM_PCA(all_mats, ...
        all_behav, thresh, v_alpha, lambda, k, seed)
    %ridgeCPM Connectome-based predictive modeling using univariate
    %feature selection and ridge regression with additional PCA to combine
    %measures
    %
    %   [q_s, q_s_fold, r_pearson, r_rank, y, new_behav, pca_coeff_total, coef_total, coef0_total, lambda_total] = ridgeCPM_PCA(all_mats, all_behav, thresh, v_alpha, lambda, k, seed)
    %
    %   Input:      all_mats,           connectome of all the subjects and tasks
    %                                   [rxegions x regions x subjects x tasks]
    %
    %               all_behav,          behavior of all the subjects
    %                                   [subjects x measures]
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
    %                                   validation or the index, default is 10
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
    %               new_behav,          the new behavior PCA generated in each
    %                                   fold
    %
    %               pca_coeff_total     PCA coefficients for each fold to
    %                                   combine different measures
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
    
    %% initialization
    if nargin < 3
        error('not enough arguments, please check the help')
    end
    
    if ~exist('k', 'var')
        k = 10;
    elseif numel(k) > 1
        if numel(k) ~= size(all_mats, 3)
            error('invalid cross validation index')
        else
            indices = k;
            k = numel(unique(indices));
        end
    end
    
    if numel(k) > 1
        if numel(k) ~= size(all_mats, 3)
            error('invalid cross validation index')
        else
            indices = k;
            k = numel(unique(indices));
        end
    end
    
    if ~exist('seed', 'var') || length(v_alpha) ~= 1 
        seed = 665;
    end
    
    if ~exist('v_alpha', 'var') || length(v_alpha) ~= 1 
        v_alpha = 1e-6;
    end
    
    num_sub_total = size(all_mats, 3);
    num_node = size(all_mats, 1);
    num_task = size(all_mats, 4);
    num_behav = size(all_behav, 2);

    is_sym = issymmetric(all_mats(:, :, 1, 1));
    if is_sym
        num_edge = num_node * (num_node - 1) / 2;
    else
        num_edge = num_node * num_node;
    end

    coef_total = zeros(num_edge*num_task, k); %store all the coefficients
    coef0_total = zeros(1, k); % store all the intercept
    lambda_total = zeros(1, k); % store all the lambda
    q_s_fold = zeros(1, k);
    pca_coef_total = zeros(k, num_behav);
    new_behav = zeros(num_sub_total, 1);

    % set cross validation index
    rng(seed, 'twister');
    if ~exist('indices', 'var')
        indices = crossvalind('Kfold', num_sub_total, k);
    end
    
    %% convert connectivity to edge matrix (could made easier by squareform)
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
    

    %% main
    y = zeros(num_sub_total, 1);   
    for i_fold = 1 : k
        tic
        fprintf('%dth fold\n', i_fold);
        
        test_idx = (indices==i_fold);
        train_idx = (indices~=i_fold);
        train_mats = all_edges(:, train_idx);
        train_behav = all_behav(train_idx, :);
        test_mats = all_edges(:, test_idx);
        test_behav = all_behav(test_idx, :);
        
        [pca_coef_total(i_fold, :)] = pca(train_behav, 'NumComponents', 1);
        new_train_behav = train_behav * pca_coef_total(i_fold, :)';
        new_test_behav = test_behav * pca_coef_total(i_fold, :)';
        new_behav(test_idx) = new_test_behav;

        % first step univariate edge selection
        [~, edge_p] = corr(train_mats', new_train_behav);
        edges_1 = find(edge_p < thresh);
        disp(numel(edges_1))
        
        % build model on TRAIN subs
        if ~exist('lambda', 'var') 
            [fit_coef, fit_info] = lasso(train_mats(edges_1, :)', new_train_behav, 'Alpha',v_alpha, 'CV', 10);
            idxLambda1SE = fit_info.Index1SE;
            coef = fit_coef(:,idxLambda1SE);
            coef0 = fit_info.Intercept(idxLambda1SE);
            lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
        elseif numel(lambda) > 1
            [fit_coef, fit_info] = lasso(train_mats(edges_1, :)', new_train_behav, 'Alpha',v_alpha, 'Lambda', lambda);
            idxLambda = find(fit_info.MSE==min(fit_info.MSE));
            coef = fit_coef(:, idxLambda);
            coef0 = fit_info.Intercept(idxLambda);
            lambda_total(i_fold) = fit_info.Lambda(idxLambda);
        else
            [coef, fit_info] = lasso(train_mats(edges_1, :)', new_train_behav, 'Alpha',v_alpha, 'Lambda', lambda);
            coef0 = fit_info.Intercept;
        end

        % run model on TEST sub with the best lambda parameter
        
        y(test_idx) = test_mats(edges_1, :)'*coef+coef0;
        
        coef_total(edges_1, i_fold) = coef;
        coef0_total(:, i_fold) = coef0;
        
        mse = sum((y(test_idx) - new_test_behav).^2) / sum(test_idx);
        q_s_fold(i_fold) = 1 - mse / var(new_test_behav, 1); 
        toc
    end
    
    % compare predicted and observed behaviors
    [r_pearson, ~] = corr(y, new_behav);
    [r_rank, ~] = corr(y, new_behav, 'type', 'spearman');
    mse = sum((y - new_behav).^2) / num_sub_total;
    q_s = 1 - mse / var(new_behav, 1);
end

