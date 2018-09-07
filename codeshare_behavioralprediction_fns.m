% Copyright 2015 Xilin Shen and Emily Finn 

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,
% Papademetris X & Constable RT. (2015). Functional connectome
% fingerprinting: Identifying individuals using patterns of brain
% connectivity. Nature Neuroscience 18, 1664-1671.

% This code provides a framework for implementing functional
% connectivity-based behavioral prediction in a leave-one-subject-out
% cross-validation scheme, as described in Finn, Shen et al 2015 (see above
% for full reference). The first input ('all_mats') is a pre-calculated
% MxMxN matrix containing all individual-subject connectivity matrices,
% where M = number of nodes in the chosen brain atlas and N = number of
% subjects. Each element (i,j,k) in these matrices represents the
% correlation between the BOLD timecourses of nodes i and j in subject k
% during a single fMRI session. The second input ('all_behav') is the
% Nx1 vector of scores for the behavior of interest for all subjects.

% As in the reference paper, the predictive power of the model is assessed
% via correlation between predicted and observed scores across all
% subjects. Note that this assumes normal or near-normal distributions for
% both vectors, and does not assess absolute accuracy of predictions (only
% relative accuracy within the sample). It is recommended to explore
% additional/alternative metrics for assessing predictive power, such as
% prediction error sum of squares or prediction r^2.


%clear;
%clc;

% ------------ INPUTS -------------------

res_struct=select_randsubs(ipmats_res, double(pmatvals'), 500, 1, 0.01);

%a = cpm_cv(ipmats_res(:,1:100), double(pmatvals(1:100)'), 10, 0.01);


% % compare predicted and observed scores
% 
% [R_pos, P_pos] = corr(behav_pred_pos,all_behav);
% [R_neg, P_neg] = corr(behav_pred_neg,all_behav);
% 
% figure(1); plot(behav_pred_pos,all_behav,'r.'); 
% figure(2); plot(behav_pred_neg,all_behav,'b.'); 

function res_struct=select_randsubs(ipmats, behav, numsubs, numiters, thresh)

    res_struct=struct();

    for iter = 1:numiters
        
        fprintf('\n Performing iter # %6.0f of %6.0f \n',iter,numiters);
        randinds=randperm(numsubs);
        randipmats=ipmats(:,randinds);
        randbehav=behav(randinds);

        % LOOCV
        res_struct.loo(iter,:) = cpm_cv(randipmats, randbehav, numsubs, thresh);
        % Split half
        res_struct.k2(iter,:) = cpm_cv(randipmats, randbehav, 2, thresh);
        % K = 5
        res_struct.k5(iter,:) = cpm_cv(randipmats, randbehav, 5, thresh);
        % K = 10
        res_struct.k10(iter,:) = cpm_cv(randipmats, randbehav, 10, thresh);
        
        
    end

end


function [fit_pos,fit_neg, pos_mask, neg_mask] = train_cpm(train_vcts,train_behav,thresh)

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


function [Rpos,Rneg,Ppos,Pneg] = cpm_cv(ipmats, behav, kfolds, thresh)
    
    nsubs=size(ipmats,2);
    randinds=randperm(nsubs);
    ksample=floor(nsubs/kfolds);
   
    
    behav_pred_pos=zeros(kfolds,ksample);
    behav_pred_neg=zeros(kfolds,ksample);
    test_behav_gather=zeros(kfolds,ksample);
    
    for leftout = 1:kfolds
        fprintf('\n Performing fold # %6.0f of %6.0f \n',leftout,kfolds);

        if kfolds == nsubs
            testinds=randinds(leftout);
            traininds=setdiff(randinds,testinds);            
        else
            si=1+((leftout-1)*ksample);        
            fi=si+ksample-1;

            testinds=randinds(si:fi);
            traininds=setdiff(randinds,testinds);
        
        end
        
        % leave out subject from matrices and behavior

        train_vcts = ipmats(:,traininds);
        train_behav = behav(traininds);
        
        [fit_pos,fit_neg,pos_mask,neg_mask] = train_cpm(train_vcts, train_behav,thresh);

        test_vcts = ipmats(:,testinds);
        test_behav = behav(testinds);
        test_behav_gather(leftout,:)=test_behav;
        
        test_sumpos = sum(test_vcts.*pos_mask)/2;
        test_sumneg = sum(test_vcts.*neg_mask)/2;

        behav_pred_pos(leftout,:) = fit_pos(1)*test_sumpos + fit_pos(2);
        behav_pred_neg(leftout,:) = fit_neg(1)*test_sumneg + fit_neg(2);
        
        
    end

    behav_pred_pos=reshape(behav_pred_pos,[],1);
    behav_pred_neg=reshape(behav_pred_neg,[],1);
    test_behav_gather=reshape(test_behav_gather,[],1);
    
    [Rpos,Ppos]=corr(test_behav_gather,behav_pred_pos);
    [Rneg,Pneg]=corr(test_behav_gather,behav_pred_neg);
    
    
end