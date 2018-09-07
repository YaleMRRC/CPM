function [fit_pos,fit_neg, pos_mask, neg_mask] = train_cpm(train_vcts,train_behav,thresh)

    % Accepts train_vcts, the feature vector for each subject, in a matrix
    % of size numfeats x nsubs
    % train_behav, behaviorual trait of interest for each subject of size
    % nsubs x 1
    % thresh, the p value to use to threshold correlations


    no_sub = size(train_vcts,2);
    no_node = sqrt(size(train_vcts,1));

    % correlate all edges with behavior
    
    [r_mat,p_mat] = corr(train_vcts',train_behav);
    
    r_mat = reshape(r_mat,no_node,no_node);
    p_mat = reshape(p_mat,no_node,no_node);
    
    % set threshold and define masks
    
    pos_mask = zeros(no_node,no_node);
    neg_mask = zeros(no_node,no_node);
    
    pos_edges = find(r_mat > 0 & p_mat < thresh);
    neg_edges = find(r_mat < 0 & p_mat < thresh);
    
    pos_mask(pos_edges) = 1;
    neg_mask(neg_edges) = 1;
    
    pos_mask=reshape(pos_mask,[],1);
    neg_mask=reshape(neg_mask,[],1);
    
    % get sum of all edges in TRAIN subs (divide by 2 to control for the
    % fact that matrices are symmetric)
    
    train_sumpos = zeros(no_sub,1);
    train_sumneg = zeros(no_sub,1);
    
    
    for ss = 1:size(train_sumpos)
        train_sumpos(ss) = sum(train_vcts(:,ss).*pos_mask)/2;
        train_sumneg(ss) = sum(train_vcts(:,ss).*neg_mask)/2;
    end
    
    % build model on TRAIN subs
    fit_pos = polyfit(train_sumpos, train_behav,1);
    fit_neg = polyfit(train_sumneg, train_behav,1);
    
end
