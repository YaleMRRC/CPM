classdef phenotype < handle
    % phenotypes = {'wms','wais'}
    properties
        name;
        all_behav; % N* T labels for each phenotype
        pca_behav; % K* T
    end
    methods
        function this = phenotype(name,all_behav)
            if nargin > 1
                this.all_behav= all_behav;
                this.name = name;
                this.pca_behav = zeros(size(all_behav,1),size(all_behav,2));
            end
        end
    end
end
