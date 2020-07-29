classdef predictory < handle
    % phenotypes = {'wms','wais'}
    properties
        group; % all subjects
        phenotype; % phenotypes
        thresh; % p-value threshold used for univariate
        %                                   feature selection
        v_alpha; % value of the alpha parameter in elastic
        %                                   net, default is 1e-6 which makes the
        %                                   regression method to be ridge
        %                                   regression, v_alpha=1 makes it lasso.
        lambda; % value of the lambda, if not provided, cross-validation will be used
        k; % (optional) number of folds in k-fold cross validation, default is 10
        seed;% (optional),     random seed, default is 665
        q_s;% cross-validated R^2 between predicted value with ground truth (all_behav)
        r_pearson;% WRONG! direct pearson correlation between predicted value with ground
        %                                   truth, only kept here for comparison will be removed afterwards
        %
        r_rank; % cross-validated spearman correlation
        %                                   between predicted value with ground truth (all_behav)
        Y;% predicted value for all the subjects
        coef_total; % regression coefficients of all the edges in all the k folds
        coef0_total;% regression intercept in all the k folds
        lambda_total; % penalty parameter chosen at each iteration
        num_sub_total;
        num_node;
        num_task;
        num_edge;
        all_edges;
        all_pos_edges; % iter*edges TODO
        all_neg_edges; % iter*edges TODO

        output;
        mse;
        control;
    end
    methods
        function this = predictory(group,options)
            this.k = 2;% default folds
            this.v_alpha = 1e-6;
            this.seed = 665;
            this.thresh=0.05;
            this.control = options.control;
            
            this.group=group;
            this.phenotype=options.phenotype;
            if length(fieldnames(options))< 3
                error('not enough arguments, please check the help')
            end
            
            if isfield(options,'k')
                this.k=options.k;
            end
            
            if isfield(options,'thresh')
                this.thresh=options.thresh;
            end
            
            if isfield(options,'seed')
                this.seed = options.seed; %Change this to put in a randomly generated number
            end
            
            if isfield(options,'v_alpha')
                this.v_alpha = options.v_valpha;
            end
            
            if isfield(options,'lambda')
                this.lambda = options.lambda;
            end
            
            this.num_node = group.num_node;
            this.num_task = group.num_task;
            this.num_edge = group.num_edge;
            this.num_sub_total = group.group_size;
            
            this.coef_total = zeros(this.num_edge*this.num_task, this.k); %store all the coefficients
            this.coef0_total = zeros(1, this.k); % store all the intercept
            this.lambda_total = zeros(1, this.k); % store all the lambda
            this.all_edges = this.group.all_edges;
            this.Y = zeros(group.group_size,1);
        end
    end
    methods (Abstract)
        run(this);
        evaluate(this);
    end
    methods
        function [train_behav,test_behav] = kfolds_pca(this,train,test,i_fold)
            %Define Mu, PCA coeff for the train set
            [coeff,score,latent,tsquared,explained,mu] = pca(train.y);
            this.output(i_fold).coeff = coeff;
            this.output(i_fold).score=score;
            this.output(i_fold).latent=latent;
            this.output(i_fold).tsquared=tsquared;
            this.output(i_fold).explained=explained;
            this.output(i_fold).mu = mu;
            train_pca = this.output(i_fold).score;
            %Apply Mu, PCA coeff on the test set (to compute "actual" PCA score)
            test_pca = (test.y-this.output(i_fold).mu)*this.output(i_fold).coeff;
            %Save these back into the original all_behav indicies in case we want to compute PCA stability across i_folds
            this.output(i_fold).allscores = zeros(size(this.phenotype.all_behav));
            this.output(i_fold).allscores(train.indx,:) = train_pca;
            this.output(i_fold).allscores(test.indx,:) = test_pca;
            train_behav = train_pca;%(:,1) Select only the first component score, which expains the greatest amt of variance.
            test_behav = test_pca;%(:,1);
        end
    end
end
