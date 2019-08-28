classdef group
    properties
        subjects; % N* 268*268
%         phenotypes;
        group_size;% N
        num_node;
        num_task;
        num_edge;
        all_edges;
        gender;
        issym;
        k_fold; % number of folds
    end
    methods
        function this = group(subjects)
            this.subjects = subjects;
            this.group_size = size(subjects,2);
            this.num_node = subjects(1).num_node;
            this.num_edge = subjects(1).num_edge;
            this.num_task = subjects(1).num_task;
            this.issym = subjects(1).issym;
            this.all_edges = zeros(this.num_edge,this.group_size,this.num_task);
            this.gender = zeros(this.group_size,1);
            for i=1:this.group_size
                this.all_edges(:,i,:) = this.subjects(i).all_edges;
            end
        end
    end
end