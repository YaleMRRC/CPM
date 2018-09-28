function [q_s, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total,output_pca,pcascores_test_all] = ridgeCPM(all_mats, all_behav, thresh, seed, v_alpha, lambda, k)
    %ridgeCPM Connectome-based predictive modeling using univariate
    %feature selection and ridge regression
    %
    %   [q_s, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total] = ridgeCPM(all_mats, all_behav, thresh, v_alpha, lambda, k, seed)
    %
    %   Input:      all_mats,           connectome of all the subjects and tasks
    %                                   [regions x regions x subjects x tasks]
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
    
    %% initialization
    if nargin < 3
        error('not enough arguments, please check the help')
    end
    
    if ~exist('k', 'var')
        k = 10;
    end
    
    if ~exist('seed', 'var')
        seed = 665; %Change this to put in a randomly generated number
    end
    
    if ~exist('v_alpha', 'var') || length(v_alpha) ~= 1 
        v_alpha = 1e-6;
    end
    
    num_sub_total = size(all_mats, 3);
    num_node = size(all_mats, 1);
    num_task = size(all_mats, 4);

    is_sym = issymmetric(all_mats(:, :, 1, 1));
    if is_sym
        num_edge = num_node * (num_node - 1) / 2;
    else
        num_edge = num_node * num_node;
    end

    coef_total = zeros(num_edge*num_task, k); %store all the coefficients
    coef0_total = zeros(1, k); % store all the intercept
    lambda_total = zeros(1, k); % store all the lambda

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
    rng(seed);
    indices = crossvalind('Kfold', num_sub_total, k);
    pcascores_test_all = zeros(size(all_behav));
    for i_fold = 1 : k
        fprintf('%dth fold\n', i_fold);
        
        test_idx = (indices==i_fold);
        train_mats = all_edges;
        train_mats(:, test_idx) = [];
        %Run PCA on this fold's data
        [output_pca(i_fold).coeff,output_pca(i_fold).score,output_pca(i_fold).latent,output_pca(i_fold).tsquared,output_pca(i_fold).explained,output_pca(i_fold).mu,test_pca,train_behav] = kfolds_PCA(all_behav,indices,i_fold);
        pcascores_test_all(test_idx,:)=test_pca;
        
        % first step univariate edge selection
        [~, edge_p] = corr(train_mats', train_behav);
        edges_1 = find(edge_p < thresh);
        
        
        % build model on TRAIN subs
        if ~exist('lambda', 'var')  || length(lambda) ~= 1 
            [fit_coef, fit_info] = lasso(train_mats(edges_1, :)', train_behav, 'Alpha',v_alpha, 'CV', 10);
            idxLambda1SE = fit_info.Index1SE;
            coef = fit_coef(:,idxLambda1SE);
            coef0 = fit_info.Intercept(idxLambda1SE);
            lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
        else
            [coef, fit_info] = lasso(train_mats(edges_1, :)', train_behav, 'Alpha',v_alpha, 'Lambda', lambda);
            coef0 = fit_info.Intercept;
        end

        % run model on TEST sub with the best lambda parameter
        test_mats = all_edges(:, test_idx);
        y(test_idx) = test_mats(edges_1, :)'*coef+coef0;
        
        coef_total(edges_1, i_fold) = coef;
        coef0_total(:, i_fold) = coef0;
        
    end
    
    % compare predicted (y) and observed (pcascores_test_all) behaviors
    [r_pearson, ~] = corr(y, pcascores_test_all);
    [r_rank, ~] = corr(y, pcascores_test_all, 'type', 'spearman');
    mse = sum((y - pcascores_test_all).^2) / num_sub_total;
    q_s = 1 - mse / var(pcascores_test_all, 1);
end

function [coeff,score,latent,tsquared,explained,mu,test_pca,train_behav] = kfolds_PCA(all_behav,indices,i_fold)
    %Define the train and test sets for i_fold
    test_idx = (indices==i_fold);
    train_behav = all_behav;
    train_behav(test_idx,:) = [];
    test_behav = all_behav(test_idx,:);
    
    %Define Mu, PCA coeff for the train set
    [coeff,score,latent,tsquared,explained,mu] = pca(train_behav);
    train_pca = score;
    
    %Apply Mu, PCA coeff on the test set (to compute "actual" PCA score)
    test_pca = (test_behav-mu)*coeff;
    test_pca = test_pca(:,1);
    %Save these back into the original all_behav indicies in case we want to compute PCA stability across i_folds 
    %pcascores_test = zeros(size(all_behav));
    %pcascores_test(indices==i_fold,:) = test_pca;
     
    train_behav = train_pca(:,1); %Select only the first component score, which expains the greatest amt of variance.
end