classdef subject < handle
    properties
        id;
        dataset;
        all_edges; % K * 268*268
        Matrix; % 268*268 matrix based on current task
        gender; % male: 1, female: 0
        age; %
        dim; % 268*268
        issym;
        num_node;
        num_task;
        num_edge;
    end
    methods
        function this = subject(x,id,dataset,mask)
            if nargin > 1
                this.id = id;
                this.dataset = dataset;
                this.num_node = size(x, 1);
                this.num_task = size(x, 4);
                is_sym = issymmetric(x(:, :, 1, 1));
                if is_sym
                    this.num_edge = this.num_node * (this.num_node - 1) / 2;
                else
                    this.num_edge = this.num_node * this.num_node;
                end
                % apply mask if it has been set
                x = squeeze(x);
                this.all_edges = zeros(this.num_edge,this.num_task);
                for j_task = 1 : this.num_task
                    if is_sym
                        this.all_edges(:, j_task) = squareform(tril(x(:, :, j_task), -1));
                    else
                        this.all_edges(:, j_task) = reshape(x(:, :, j_task), [], 1);
                    end
                end
           
                this.issym = issymmetric(x(:, :, 1, 1));
                if(mask=="abi")
                    missing_nodes = [60, 100, 108, 109, 112, 115,116, 118, 129, 189, ...
                        202, 239, 240, 242, 243, 249, 250, 266];
                    f1 = load('../data/abby_all_task_edges_pos_p0.001.mat');
                    f2 = load('../data/abby_all_task_edges_neg_p0.001.mat');
                    mask = squeeze(f1.task_pos_edge(2,:,:)+f2.task_neg_edge(2,:,:)); % wm mask
                    for i=missing_nodes
                        b1 = zeros(1,size(mask,1));
                        b2 = zeros(size(mask,1)+1,1);
                        mask = [mask(1:i,:); b1; mask(i+1:end,:)];
                        mask = [mask(:,1:i) b2 mask(:,i+1:end)];
                    end % now we have 268*268
                elseif(mask=="nbs" && contains(this.dataset,"hcp"))
                    mask= csvread('../data.515/nbs.100.hcp.all_tasks.csv');
                    this.all_edges = zeros(this.num_edge,this.num_task);
                    for j_task = 1 : this.num_task
                        mat = x(:, :, j_task).*mask;
                        this.all_edges(:, j_task) = squareform(mat);
                    end
                end
            end
        end
    end
end