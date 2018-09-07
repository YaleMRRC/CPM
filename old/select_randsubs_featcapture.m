function [res_struct,behav_struct]=select_randsubs_featcapture(ipmats, behav, numsubs, numiters, thresh, ipmats_ex,external_compare,behav_ex, normalize)

    res_struct=struct();

    behav_struct=struct();
    
    behav_popvar=mean((behav-mean(behav)).^2);

    
    for iter = 1:numiters
        
        fprintf('\n Performing iter # %6.0f of %6.0f \n',iter,numiters);
        totalsubs=size(ipmats,2);
        randinds=randperm(totalsubs);
        randinds=randinds(1:numsubs);
        randipmats=ipmats(:,randinds);
        randbehav=behav(randinds);
        res_struct.subject_inds(iter,:)=randinds;
        
        
        if normalize
            behav_mean=mean(randbehav);
            behav_var=var(randbehav);
            randbehav=(randbehav-behav_mean)/behav_var;
            
            res_struct.behav_mean(iter)=behav_mean;
            res_struct.behav_var(iter)=behav_var;
            
            if external_compare
                behav_ex=(behav_ex-mean(behav_ex))/var(behav_ex);
            end        
        end

        % LOOCV
        [res_struct.loo(iter,1),res_struct.loo(iter,2),res_struct.loo(iter,3),res_struct.loo(iter,4),res_struct.loo(iter,5),res_struct.loo(iter,6), ...
            behav_struct.loo.testbehav(iter,:), behav_struct.loo.predbehavpos(iter,:), behav_struct.loo.predbehavneg(iter,:), ...
            behav_struct.loo.pos_feats(iter,:,:), behav_struct.loo.neg_feats(iter,:,:)] ...
            = cpm_cv(randipmats, randbehav, numsubs, thresh, behav_popvar);
        
        % Split half
        [res_struct.k2(iter,1),res_struct.k2(iter,2),res_struct.k2(iter,3),res_struct.k2(iter,4),res_struct.k2(iter,5),res_struct.k2(iter,6), ...
            behav_struct.k2.testbehav(iter,:), behav_struct.k2.predbehavpos(iter,:), behav_struct.k2.predbehavneg(iter,:), ...
            behav_struct.k2.pos_feats(iter,:,:), behav_struct.k2.neg_feats(iter,:,:)] ... 
            = cpm_cv(randipmats, randbehav, 2, thresh, behav_popvar);
        
        
        % K = 5
        [res_struct.k5(iter,1),res_struct.k5(iter,2),res_struct.k5(iter,3),res_struct.k5(iter,4),res_struct.k5(iter,5),res_struct.k5(iter,6), ...
            behav_struct.k5.testbehav(iter,:), behav_struct.k5.predbehavpos(iter,:), behav_struct.k5.predbehavneg(iter,:), ...
            behav_struct.k5.pos_feats(iter,:,:), behav_struct.k5.neg_feats(iter,:,:)] ... 
            = cpm_cv(randipmats, randbehav, 5, thresh, behav_popvar);
        
        % K = 10
        [res_struct.k10(iter,1),res_struct.k10(iter,2),res_struct.k10(iter,3),res_struct.k10(iter,4),res_struct.k10(iter,5),res_struct.k10(iter,6), ...
            behav_struct.k10.testbehav(iter,:), behav_struct.k10.predbehavpos(iter,:), behav_struct.k10.predbehavneg(iter,:), ...
            behav_struct.k10.pos_feats(iter,:,:), behav_struct.k10.neg_feats(iter,:,:)] ... 
             = cpm_cv(randipmats, randbehav, 10, thresh, behav_popvar);
        
        if external_compare
            %# External Val
            [fit_pos,fit_neg, pos_mask, neg_mask] = train_cpm(randipmats,randbehav,thresh);


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
            behav_struct.external.pos_feats(iter,:)=pos_mask;
            behav_struct.external.neg_feats(iter,:)=neg_mask;

        end
    end
end
