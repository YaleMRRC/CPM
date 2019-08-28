classdef explanatory
    %EXPLANATORY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        group;
        phenotypes;
        diagnosis;
        thresh;
        coef_total; % regression coefficients of all the edges in all the k folds
        coef0_total;% regression intercept in all the k folds
        lambda_total; % penalty parameter chosen at each iteration
        num_sub_total;
        seed;
        lambda;
        v_alpha;
        num_node;
        num_task;
        num_edge;
        all_edges;
        output;
    end
    
    methods
        function this = explanatory(group,options)
            %EXPLANATORY Construct an instance of this class
            %   Detailed explanation goes here
            this.phenotypes = options.phenotype;
            this.diagnosis = options.diagnosis;
            this.group = group;
            
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
            
            this.lambda_total = 0; % store all the lambda
            this.all_edges = this.group.all_edges;
        end
        function [pID,pN] = my_fdr(this,p,q)
            % FORMAT [pID,pN] = FDR(p,q)
            %
            % p   - vector of p-values
            % q   - False Discovery Rate level
            %
            % pID - p-value threshold based on independence or positive dependence
            % pN  - Nonparametric p-value threshold
            %______________________________________________________________________________
            % $Id: FDR.m,v 2.1 2010/08/05 14:34:19 nichols Exp $
            p = p(isfinite(p));  % Toss NaN's
            p = sort(p(:));
            V = length(p);
            I = (1:V)';
            
            cVID = 1;
            cVN = sum(1./(1:V));
            
            pID = p(max(find(p<=I/V*q/cVID)));
            if isempty(pID), pID=0; end
            pN = p(max(find(p<=I/V*q/cVN)));
            if isempty(pN), pN=0; end
        end
        
        function [mask,biggestArea] = my_largestCC(this,a)
            %Brett's shot at Doug's "Biggest 4-connected" puzzler
            %08/18/08
            lindices = tril(ones(268));
            lindices(find(eye(268)))=0;
            lindices = find(lindices);
            
            network = zeros(268,268);
            network(lindices) = a;
            network = network'+network;
            
            biggestArea = 0;
            L = bwlabel(network==1,4);
            stats = regionprops(L,'area');
            [tmp,tmpbig] = max([stats.Area]);
            if tmp > biggestArea
                biggestArea = tmp;
                mask = L == tmpbig;
            end
        end
    end
end
