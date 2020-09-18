function [q_s, r_pearson, r_rank, y, all_edge_weight, mask] = cCPM(all_mats, ...
        all_behav, thresh, k, seed)
    %CCA-CPM Connectome-based predictive modeling using univariate feature selection 
    %
    %   [q_s, r_pearson, r_rank, y, all_edge_weight, mask] = cCPM(all_mats, all_behav, 0.1)
    %   [q_s, r_pearson, r_rank, y, all_edge_weight, mask] = cCPM(all_mats, all_behav, 0.1, 10, 213)
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
    %               k(optional),        number of folds in k-fold cross
    %                                   validation, default is 10
    %
    %               seed(optional),     random seed, default is 665
    %
    %   Output:     q_s,                cross-validated R^2 between predicted
    %                                   value with ground truth (all_behav)
    %
    %               r_pearson,          cross-validated pearson correlation 
    %                                   between predicted value with ground 
    %                                   truth (all_behav)
    %
    %               r_rank,             cross-validated spearman correlation
    %                                   between predicted value with ground 
    %                                   truth (all_behav)
    %
    %               y,                  predicted value for all the subjects
    %  
    %               all_edge_weight,    CCA projection weight of each edge
    % 
    %               mask,               network selected
    %
    %   Reference: Siyuan Gao, ISBI 2018
    %   Siyuan Gao, Yale University, 2018-2019
    
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
    
    num_sub_total = size(all_mats, 3);
    num_node = size(all_mats, 1);
    num_task = size(all_mats, 4);
    is_sym = issymmetric(all_mats(:, :, 1, 1));
    if is_sym
        num_edge = num_node * (num_node - 1) / 2;
    else
        num_edge = num_node * num_node;
    end
    
    mask = zeros(num_node, num_node, k); %store all the edges selected
    
    %% convert connectivity to edge matrix
    all_edges = zeros(num_edge, num_sub_total, num_task);
    for i_sub = 1 : num_sub_total
        for j_task = 1 : num_task
            if is_sym
                all_edges(:, i_sub, j_task) = squareform(tril(all_mats(:, :, ...
                    i_sub, j_task), -1));
            else
                all_edges(:, i_sub, j_task) = reshape(all_mats(:, :, i_sub, ...
                    j_task), [], 1);
            end
        end
    end

    
    %% main
    
    y = zeros(num_sub_total, 1);
    all_edge_weight = zeros(k, num_edge, num_task);
    
    rng(seed);
    indices = crossvalind('Kfold', num_sub_total, k);
    
    for i_fold = 1 : k
        fprintf('%dth fold\n', i_fold);
        
        test_idx = (indices==i_fold);
        train_mats = all_edges;
        train_mats(:, test_idx, :) = [];
        train_behav = all_behav;
        train_behav(test_idx, :) = [];
        num_train = size(train_behav, 1);

        % weight calculation
        edge_p = zeros(num_edge, 1);
        weight = zeros(num_edge, num_task);
        for j_edge = 1 : num_edge
            %disp(j/num_edge)
            edge_temp = squeeze(train_mats(j_edge, :, :));
            if rank(edge_temp)==num_task
                [A, B, ~, ~, ~, stats] = canoncorr(edge_temp, train_behav);
                if B > 0
                    weight(j_edge, :) = A;
                elseif B<0
                    weight(j_edge, :) = -A;
                end
                edge_p(j_edge) = stats.pF;
            end
        end

        mean_feature = squeeze(mean(train_mats, 2));
        
        % define masks
        edge_idx = edge_p < thresh;
        disp(sum(edge_idx))
        
        all_edge_weight(i_fold, edge_idx, :) = weight(edge_idx, :);
        % get sum of all edges in TRAIN subs (divide by 2 to control for the
        % fact that matrices are symmetric)
        train_sum = zeros(num_train,1);

        for j_train = 1:num_train
            train_mats_temp = sum((squeeze(train_mats(:, j_train, :)) - ...
                mean_feature) .* weight, 2);
            train_sum(j_train) = (sum(sum(train_mats_temp .* edge_idx))) / 2;
        end

        % build model on TRAIN subs
        fit_info = polyfit(train_sum, train_behav,1);

        % run model on TEST sub
        test_mats = all_edges(:, test_idx, :);
        num_test = size(test_mats, 2);
        test_sum = zeros(num_test, 1);
        
        for j_test = 1 : num_test
             test_mats_temp = sum((squeeze(test_mats(:, j_test, :)) - ...
                 mean_feature) .* weight, 2);
             test_sum(j_test) = (sum(sum(test_mats_temp .* edge_idx))) / 2;
        end
         
        y(test_idx) = fit_info(1).*test_sum + fit_info(2);
        
        % save out the selected edges in this fold
        mask(:, :, i_fold) = squareform(edge_idx);
        
    end
    
    % compare predicted and observed behaviors
    [r_pearson, ~] = corr(y, all_behav);
    [r_rank, ~] = corr(y, all_behav, 'type', 'spearman');
    mse = sum((y - all_behav).^2) / num_sub_total;
    q_s = 1 - mse / var(all_behav, 1);
end
