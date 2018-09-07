function [res_struct,pred_behav_struct]=nsubs_iter(ipmats,behav,numiters,thresh,normalize)

    res_struct=struct();

    for iter = 1:numiters
    
        fprintf('\n Performing iter # %6.0f of %6.0f \n',iter,numiters);

        
        nsubs=size(ipmats,2);
        randinds=randperm(nsubs);
        testinds=randinds(1:400);
        testmats=ipmats(:,testinds);
        testbehav=behav(testinds);

        remaining_rand=randinds(401:end);
        
        if normalize
            testbehav_mean=mean(testbehav);
            testbehav_var=var(testbehav);
            testbehav=(testbehav-testbehav_mean)/testbehav_var;
            pred_behav_struct.testbehavmean(iter)=testbehav_mean;
            pred_behav_struct.testbehavvar(iter)=testbehav_var;
        end
        
        pred_behav_struct.testbehav(iter,:)=testbehav;

        
        for trainsubs = 25:25:400
            fprintf('\n Training on %6.0f Subs \n',trainsubs);

            traininds=remaining_rand(1:trainsubs);
            trainmats=ipmats(:,traininds);
            trainbehav=behav(traininds);
            
            if normalize
                trainbehav_mean=mean(trainbehav);
                trainbehav_var=var(trainbehav);
                trainbehav=(trainbehav-trainbehav_mean)/trainbehav_var;
                pred_behav_struct.(['train' num2str(trainsubs)]).trainbehavmean(iter)=trainbehav_mean;
                pred_behav_struct.(['train' num2str(trainsubs)]).trainbehavvar(iter)=trainbehav_var;
            end
            
            
            [fit_pos,fit_neg,pos_mask,neg_mask] = train_cpm(trainmats, trainbehav,thresh);

            test_sumpos = sum(testmats.*repmat(pos_mask,1,400))/2;
            test_sumneg = sum(testmats.*repmat(neg_mask,1,400))/2;

            behav_pred_pos = fit_pos(1)*test_sumpos + fit_pos(2);
            behav_pred_neg = fit_neg(1)*test_sumneg + fit_neg(2);

            [Rpos,Ppos]=corr(testbehav,behav_pred_pos');
            [Rneg,Pneg]=corr(testbehav,behav_pred_neg');
            
            mse_pos=mean((testbehav-behav_pred_pos').^2);
            mse_neg=mean((testbehav-behav_pred_neg').^2);
            
            res_struct.(['train' num2str(trainsubs)])(iter,:)=[Rpos Rneg Ppos Pneg mse_pos mse_neg];
            
            
            pred_behav_struct.(['train' num2str(trainsubs)]).predbehavpos(iter,:)=behav_pred_pos;
            pred_behav_struct.(['train' num2str(trainsubs)]).predbehavneg(iter,:)=behav_pred_neg;
            
        end    
        
    end


end