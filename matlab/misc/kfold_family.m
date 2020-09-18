function id = kfold_family(family_id, k, seed)
    % generate k-fold idx considering familty structure
    % Inputs: family_id [n_sub x 2], first column is subject id while second
    %         column is the family id
    %        
    %         k, number of folds
    %
    %         seed, random seed used
    %
    % Outputs: id [n_sub x 1], k-fold index for each subject
    
    if ~exist('k', 'var')
        k = 10;
    end
    
    if ~exist('seed', 'var')
        seed = 213;
    end
    
    if length(unique(family_id(:,2))) < k
        error('fold number is bigger than family number, please reduce k!')
    end
    n_sub = size(family_id, 1);
    fold_size = floor(n_sub / k);
    
    id = zeros(n_sub, 1);
    
    rng(seed)
    p = randperm(n_sub);
    p_idx = 1;
    
    i_fold = 1;
    end_flag = 0;
    while i_fold <= k
        temp_size = sum(id == i_fold);
        while temp_size < fold_size && p_idx <= n_sub
            while p_idx <= n_sub && id(p(p_idx)) ~= 0
                p_idx = p_idx + 1;
            end
            if p_idx > n_sub
                end_flag = 1;
            else
                temp_family = family_id(p(p_idx), 2);
                id(family_id(:, 2)==temp_family) = i_fold;
                temp_size = sum(id == i_fold);
            end
        end
        if end_flag
            break
        end
        i_fold = i_fold + 1;
    end
end