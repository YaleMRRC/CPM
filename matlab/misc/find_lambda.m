function [opt_lambda, mse_fold] = find_lambda(all_mats, all_behav, lambda,...
        thresh, n_fold)
    lambda = sort(lambda(:),1,'descend'); % sort the lambda
    if ~exist('n_fold', 'var') || numel(n_fold) < 1
        n_fold = 2; % by-default, split half to speed up
    end
    n_lambda = numel(lambda);
    indices = crossvalind('Kfold', size(all_mats, 2), n_fold);
    mse_fold = zeros(n_fold, n_lambda);
    parfor i_lam = 1 : n_lambda
        y_fit_temp = zeros(size(all_mats, 2), 1);
        for j_fold = 1 : n_fold
%             disp(['#inner fold:', num2str(j_fold)])
            train_mats = all_mats(:, indices~=j_fold);
            train_behav = all_behav(indices~=j_fold);
            test_mats = all_mats(:, indices==j_fold);
            test_behav = all_behav(indices==j_fold);
            
            [~, edge_p] = corr(train_mats', train_behav);
            edge_idx = find(edge_p <= thresh);
            if size(train_mats,1)<300 % if few features (e.g., node activation), keep defaults
                [fit_coef, fit_info] = lasso(train_mats(edge_idx, :)',...
                    train_behav, 'Alpha',1e-9, 'Lambda', lambda(i_lam));
            else % if many features, increase RelTol and decrease numiters
                [fit_coef, fit_info] = lasso(train_mats(edge_idx, :)',...
                    train_behav, 'Alpha',1e-9, 'Lambda', lambda(i_lam),'RelTol',1e-3,'MaxIter',100);
            end
            y_fit_temp(indices==j_fold) = test_mats(edge_idx, :)'*fit_coef+...
                fit_info.Intercept;
            mse_fold(j_fold, i_lam) = mean((y_fit_temp(indices==j_fold)-...
                test_behav).^2);
        end
    end
    mean_mse = mean(mse_fold);
    se_mse = std(mse_fold) ./ sqrt(n_fold);
    min_idx = find(mean_mse==min(mean_mse));
    min_lambda = lambda(min_idx);
    min_plus1 = mean_mse(min_idx) + se_mse(min_idx);
    se_idx = find((mean_mse <= min_plus1),1,'first');
    se_lambda = lambda(se_idx);
    if se_idx == 1 % if the biggest lambda is chosen, use min_lambda
        opt_lambda = min_lambda;
    else % otherwise choose a 1se lambda that gives bigger penalty
        opt_lambda = se_lambda;
    end
end
    
            
            

    
