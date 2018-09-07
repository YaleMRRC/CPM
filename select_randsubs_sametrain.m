function [res_struct,behav_struct]=select_randsubs_sametrain(ipmats, behav, numtrain, numiters, thresh, ipmats_ex, behav_ex, normalize)

    res_struct=struct();
    
    behav_struct=struct();
    
    behav_popvar=mean((behav-mean(behav)).^2);

    for iter = 1:numiters
        
        
        
        looind=numtrain+1;
        k2ind=numtrain/0.5;
        k5ind=numtrain/0.8;
        k10ind=numtrain/0.9;
        
        
        
        fprintf('\n Performing iter # %6.0f of %6.0f \n',iter,numiters);
        totalsubs=size(ipmats,2);
        randinds=randperm(totalsubs);
        randinds=randinds(1:k2ind);
        randipmats=ipmats(:,randinds);
        randbehav=behav(randinds);

        randbehav_loo=randbehav(1:looind);
        randbehav_k2=randbehav(1:k2ind);
        randbehav_k5=randbehav(1:k5ind);
        randbehav_k10=randbehav(1:k10ind);

        
        behav_struct.loo.subject_inds(iter,:)=randinds(1:looind);
        behav_struct.k2.subject_inds(iter,:)=randinds;
        behav_struct.k5.subject_inds(iter,:)=randinds(1:k5ind);
        behav_struct.k10.subject_inds(iter,:)=randinds(1:k10ind);
        
        if normalize
                        
            % LOO
            behav_mean_loo=mean(randbehav_loo);
            behav_var_loo=var(randbehav_loo);
            randbehav_loo=(randbehav_loo-behav_mean_loo)/behav_var_loo;
            
            behav_struct.loo.behav_mean(iter)=behav_mean_loo;
            behav_struct.loo.behav_var(iter)=behav_var_loo;
            
            % K2
            behav_mean_k2=mean(randbehav_k2);
            behav_var_k2=var(randbehav_k2);
            randbehav_k2=(randbehav_k2-behav_mean_k2)/behav_var_k2;
            
            behav_struct.k2.behav_mean(iter)=behav_mean_k2;
            behav_struct.k2.behav_var(iter)=behav_var_k2;
            
            % K5
            behav_mean_k5=mean(randbehav_k5);
            behav_var_k5=var(randbehav_k5);
            randbehav_k5=(randbehav_k5-behav_mean_k5)/behav_var_k5;
            
            behav_struct.k5.behav_mean(iter)=behav_mean_k5;
            behav_struct.k5.behav_var(iter)=behav_var_k5;
            
            % k10
            behav_mean_k10=mean(randbehav_k10);
            behav_var_k10=var(randbehav_k10);
            randbehav_k10=(randbehav_k10-behav_mean_k10)/behav_var_k10;
            
            behav_struct.k10.behav_mean(iter)=behav_mean_k10;
            behav_struct.k10.behav_var(iter)=behav_var_k10;
                        
            behav_ex=(behav_ex-mean(behav_ex))/var(behav_ex);
            
        end
        
        
        % LOOCV
        [res_struct.loo(iter,1),res_struct.loo(iter,2),res_struct.loo(iter,3),res_struct.loo(iter,4),res_struct.loo(iter,5),res_struct.loo(iter,6), ...
            behav_struct.loo.testbehav(iter,:), behav_struct.loo.predbehavpos(iter,:), behav_struct.loo.predbehavneg(iter,:)] ...
            = cpm_cv(randipmats(:,1:looind), randbehav_loo, looind, thresh,behav_popvar);
        % Split half
        [res_struct.k2(iter,1),res_struct.k2(iter,2),res_struct.k2(iter,3),res_struct.k2(iter,4),res_struct.k2(iter,5),res_struct.k2(iter,6), ...
            behav_struct.k2.testbehav(iter,:), behav_struct.k2.predbehavpos(iter,:), behav_struct.k2.predbehavneg(iter,:)] ... 
            = cpm_cv(randipmats(:,1:k2ind), randbehav_k2, 2, thresh,behav_popvar);
        % K = 5
        [res_struct.k5(iter,1),res_struct.k5(iter,2),res_struct.k5(iter,3),res_struct.k5(iter,4),res_struct.k5(iter,5),res_struct.k5(iter,6), ...
            behav_struct.k5.testbehav(iter,:), behav_struct.k5.predbehavpos(iter,:), behav_struct.k5.predbehavneg(iter,:)] ...
            = cpm_cv(randipmats(:,1:k5ind), randbehav_k5, 5, thresh,behav_popvar);
        % K = 10
        [res_struct.k10(iter,1),res_struct.k10(iter,2),res_struct.k10(iter,3),res_struct.k10(iter,4),res_struct.k10(iter,5),res_struct.k10(iter,6), ...
            behav_struct.k10.testbehav(iter,:), behav_struct.k10.predbehavpos(iter,:), behav_struct.k10.predbehavneg(iter,:)] ...
             = cpm_cv(randipmats(:,1:k10ind), randbehav_k10, 10, thresh,behav_popvar);
        % External Val
        [fit_pos,fit_neg, pos_mask, neg_mask] = train_cpm(randipmats(:,1:numtrain),randbehav(1:numtrain),thresh);
        
        
        nsubs_ex=size(ipmats_ex,2);
        
        % Sum external edges
        ext_sumpos = sum(ipmats_ex.*repmat(pos_mask,1,nsubs_ex))/2;
        ext_sumneg = sum(ipmats_ex.*repmat(neg_mask,1,nsubs_ex))/2;

        % Derive predicted values for external dataset
        behav_pred_pos_ext = fit_pos(1)*ext_sumpos + fit_pos(2);
        behav_pred_neg_ext = fit_neg(1)*ext_sumneg + fit_neg(2);
        
        % Correlated predicted and real values in external datset
        [Rpos_ext,Ppos_ext]=corr(behav_ex,behav_pred_pos_ext');
        [Rneg_ext,Pneg_ext]=corr(behav_ex,behav_pred_neg_ext');
        
        % Estimate correlation based on MSE
        behav_popvar_ext=mean((behav_ex-mean(behav_ex)).^2);
        
        mse_pos=mean((behav_ex-behav_pred_pos_ext').^2);
        mse_neg=mean((behav_ex-behav_pred_neg_ext').^2);
    
        Rmsepos=sqrt(1-mse_pos/behav_popvar_ext);
        Rmseneg=sqrt(1-mse_neg/behav_popvar_ext);
        
        % Record external correlation values
        res_struct.external(iter,:) = [Rpos_ext, Rneg_ext, Ppos_ext, Pneg_ext, Rmsepos, Rmseneg];
        
        % Record test and predicted behaviour from external dataset testing
        behav_struct.external.predbehavpos(iter,:)=behav_pred_pos_ext;
        behav_struct.external.predbehavneg(iter,:)=behav_pred_neg_ext;
        behav_struct.external.testbehav(iter,:)=behav_ex;
        
    end

end
