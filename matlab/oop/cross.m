classdef cross < predictory
        properties
            group1;
            group2;
            phenotype1;
            phenotype2;
            model;
        end
    methods
        function this = cross(group1,group2,options)
            this = this@predictory(group1,options);
            this.group1=group1;
            this.group2 = group2;
            if isfield(options,'taskID1')
                if options.taskID1
                    this.group1.all_edges = squeeze(this.group1.all_edges(:,:,options.taskID1));
                else
                    xp = permute(this.group1.all_edges,[1 3 2]);
                    this.group1.all_edges = reshape(xp,[],this.group1.group_size);
                end
            end
            this.model=options.model;
            if isfield(options,'taskID2')
                if options.taskID2
                    this.group2.all_edges = squeeze(this.group2.all_edges(:,:,options.taskID2));
                else
                    xp = permute(this.group2.all_edges,[1 3 2]);
                    this.group2.all_edges = reshape(xp,[],this.group2.group_size);
                end
            end
            if(size(this.group1.all_edges,1)~=size(this.group2.all_edges,1))
                error('Number of tasks in two groups should be the same!')
            end
            if isfield(options,'phenotype1')
                this.phenotype1=options.phenotype1;
            end
            if isfield(options,'phenotype2')
                this.phenotype2=options.phenotype2;
            end

        end
        function run(this)
            if(this.model=="cpm")
                this.run_cpm();
            elseif(this.model=="rcpm")
                this.run_rcpm();
            end
        end
        function run_cpm(this)
            test.x = this.group2.all_edges;
            train.x = this.group1.all_edges;
            test.y = this.phenotype2.all_behav;
            train.y = this.phenotype1.all_behav;
            %
            % first step univariate edge selection
            [edge_corr, edge_p] = corr(train.x', train.y);
            edges_pos = (edge_p < this.thresh) & (edge_corr > 0);
            edges_neg = (edge_p < this.thresh) & (edge_corr < 0);
            
            % build model on TRAIN subs
            train_sum = (sum(train.x(edges_pos, :), 1) -sum(train.x(edges_neg, :), 1))';
            fit_train = polyfit(train_sum, train.y, 1);
            
            % run model on TEST sub
            test_sum = sum(test.x(edges_pos, :), 1) - sum(test.x(edges_neg), 1);
            this.Y = (test_sum*fit_train(1)+fit_train(2))';
        end
        function run_rcpm(this)
            all_edges = this.group2.all_edges;
            all_edges = permute(all_edges, [1, 3, 2]);
            all_edges = reshape(all_edges, [],this.group2.group_size);
            
            test.x = this.group2.all_edges;
            train.x = this.group1.all_edges;
            test.y = this.phenotype2.all_behav;
            train.y = this.phenotype1.all_behav;
            
            [~, edge_p] = corr(train.x', train.y);
            edges_1 = find(edge_p < this.thresh);

            % build model on TRAIN subs
            if ~exist('lambda', 'var')  || length(this.lambda) ~= 1
                [fit_coef, fit_info] = lasso(train.x(edges_1, :)', train.y, 'Alpha',this.v_alpha, 'CV', 10);
                idxLambda1SE = fit_info.Index1SE;
                coef = fit_coef(:,idxLambda1SE);
                coef0 = fit_info.Intercept(idxLambda1SE);
            else
                [coef, fit_info] = lasso(train.x(edges_1, :)', train.y, 'Alpha',this.alpha, 'Lambda', this.lambda);
                coef0 = fit_info.Intercept;
            end
            % run model on TEST sub with the best lambda parameter
            test.x = all_edges;
            this.Y = test.x(edges_1, :)'*coef+coef0;
        end
        function evaluate(this)
            [this.r_pearson, ~] = corr(this.Y, this.phenotype2.all_behav);
            [this.r_rank, ~] = corr(this.Y, this.phenotype2.all_behav, 'type', 'spearman');
            this.mse = sum((this.Y - this.phenotype2.all_behav).^2)/this.group2.group_size;            
            this.q_s = 1 - this.mse / var(this.Y, 1);
            fprintf('q_s=%f\n',this.q_s);
            fprintf('spearman=%f\n',this.r_rank);
            fprintf('mse=%f\n',this.mse);
            fprintf('pearson=%f\n',this.r_pearson);

        end
    end
end
