function [q_s, q_s_fold, r_pearson, r_rank, y, beta_total, intercept_total, ...
        lambda_total] = ridgeCPM(all_mats, all_behav, thresh, v_alpha,...
        lambda, k, seed)
    %ridgeCPM Connectome-based predictive modeling using univariate
    %feature selection and ridge regression
    %
    %   [q_s, q_s_fold, r_pearson, r_rank, y, coef_total, coef0_total,...
    %    lambda_total] = ridgeCPM(all_mats, all_behav, thresh, v_alpha, ...
    %    lambda, k, seed)
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
    %               beta_total,         regression coefficients of all the edges
    %                                   in all the k folds
    %
    %               intercept_total,        regression intercept in all the k folds
    %
    %               lambda_total,       penalty parameter chosen at each
    %                                   iteration
    %       
    %   Siyuan Gao, Yale University, 2018-2020
    
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
    
    if ~exist('seed', 'var')
        seed = 665;
    end
    
    if ~exist('v_alpha', 'var') || length(v_alpha) ~= 1 
        v_alpha = 1e-9;
    end
    
    % set the number of lambda to search
    numLambda = 16;
    
    % parse the dimensions of data
    num_sub_total = size(all_mats, 3);
    num_node = size(all_mats, 1);
    num_node2 = size(all_mats, 2);
    num_task = size(all_mats, 4);

    % determine the number of edges based on whether matrix is symmetric
    is_sym = issymmetric(all_mats(:, :, 1, 1));
    if is_sym
        num_edge = num_node * (num_node - 1) / 2;
    else
        num_edge = num_node * num_node2;
    end

    beta_total = zeros(num_edge*num_task, k); % store all the coefficients
    intercept_total = zeros(1, k); % store all the intercept
    lambda_total = zeros(1, k); % store all the lambda
    q_s_fold = zeros(1, k);
    
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
                all_edges(:, i_sub, j_task) = squareform(...
                    tril(all_mats(:, :, i_sub, j_task), -1));
            else
                all_edges(:, i_sub, j_task) = reshape(...
                    all_mats(:, :, i_sub, j_task), [], 1);
            end
        end
    end
    all_edges = reshape(permute(all_edges, [1, 3, 2]), [], num_sub_total);
    
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
        
        % First step univariate edge selection
        [~, edge_p] = corr(train_mats', train_behav);
        edge_idx = find(edge_p <= thresh);
        disp(['#edge: ', num2str(numel(edge_idx))])
        
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
                lambda_squence, thresh); % TIME-CONSUMING
%             disp(opt_lambda)
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
               
        beta_total(edge_idx, i_fold) = beta;
        intercept_total(:, i_fold) = intercept;
        
        mse = sum((y(test_idx) - test_behav).^2) / sum(test_idx);
        q_s_fold(i_fold) = 1 - mse / var(test_behav, 1); 
        toc
    end
    
    % compare predicted and observed behaviors
    [r_pearson, ~] = corr(y, all_behav);
    [r_rank, ~] = corr(y, all_behav, 'type', 'spearman');
    mse = sum((y - all_behav).^2) / num_sub_total;
    q_s = 1 - mse / var(all_behav, 1);
    
end

