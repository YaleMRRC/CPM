classdef rcpm < predictory
    methods
        function this = rcpm(group,options)
            this = this@predictory(group,options);
        end
        function run(this)
            
            all_edges = this.group.all_edges;
            all_edges = permute(all_edges, [1, 3, 2]);
            all_edges = reshape(all_edges, [],this.num_sub_total);
            
            rng(this.seed);
            indices = cvpartition(this.num_sub_total,'k',this.k);
            
            for i_fold = 1 : this.k
                fprintf('%dth fold\n', i_fold);
                test.indx = indices.test(i_fold);
                train.indx = indices.training(i_fold);
                test.x = this.all_edges(:,test.indx);
                train.x = this.all_edges(:,train.indx);
                test.y = this.phenotype.all_behav(indices.test(i_fold),:);
                train.y = this.phenotype.all_behav(indices.training(i_fold),:);
                
                % first step univariate edge selection
                [~, edge_p] = corr(train.x', train.y);
                edges_1 = find(edge_p < this.thresh);
                
                % build model on TRAIN subs
                if ~exist('lambda', 'var')  || length(this.lambda) ~= 1
                    [fit_coef, fit_info] = lasso(train.x(edges_1, :)', train.y, 'Alpha',this.v_alpha, 'CV', 10);
                    idxLambda1SE = fit_info.Index1SE;
                    coef = fit_coef(:,idxLambda1SE);
                    coef0 = fit_info.Intercept(idxLambda1SE);
                    this.lambda_total(i_fold) = fit_info.Lambda(idxLambda1SE);
                else
                    [coef, fit_info] = lasso(train.x(edges_1, :)', train.y, 'Alpha',this.alpha, 'Lambda', this.lambda);
                    coef0 = fit_info.Intercept;
                end
                
                % run model on TEST sub with the best lambda parameter
                test.x = all_edges(:, test.indx);
                this.y(test.indx) = test.x(edges_1, :)'*coef+coef0;
                
                this.coef_total(edges_1, i_fold) = coef;
                this.coef0_total(:, i_fold) = coef0;
            end
        end
        
        function evaluate(this)
            % compare predicted and observed behaviors
            [r_pearson, ~] = corr(this.y, this.phenotype.all_behav);
            [r_rank, ~] = corr(this.y, this.phenotype.all_behav, 'type', 'spearman');
            mse = sum((this.y - this.phenotype.all_behav).^2) / this.group.group_size;
            q_s = 1 - mse / var(this.phenotype.all_behav, 1);
            fprintf('mse=%f\n',mse);
            fprintf('r_pearson=%f\n',r_pearson);
            fprintf('r_rank=%f\n',r_rank);
            fprintf('q_s=%f\n',q_s);
        end
    end
end