classdef cpm < predictory
    methods
        function this = cpm(group,options)
            this = this@predictory(group,options);
        end
        function run(this)
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
                [edge_corr, edge_p] = corr(train.x', train.y);
                edges_pos = (edge_p < this.thresh) & (edge_corr > 0);
                edges_neg = (edge_p < this.thresh) & (edge_corr < 0);
                
                % build model on TRAIN subs
                train_sum = (sum(train.x(edges_pos, :), 1) - sum(train.x(edges_neg, :), 1))';
                fit_train = polyfit(train_sum, train.y, 1);
                
                % run model on TEST sub
                test_sum = sum(test.x(edges_pos, :), 1) - sum(test.x(edges_neg, :), 1);
                
                this.Y(test.indx) = (test_sum*fit_train(1)+fit_train(2))';
            end
        end
        function evaluate(this)
            [this.r_pearson, ~] = corr(this.Y, this.phenotype.all_behav);
            [this.r_rank, ~] = corr(this.Y, this.phenotype.all_behav, 'type', 'spearman');
            this.mse = sum((this.Y - this.phenotype.all_behav).^2) / this.num_sub_total;
            this.q_s = 1 - this.mse / var(this.Y, 1);
            fprintf('q_s=%f\n',this.q_s);
        end
    end
end
