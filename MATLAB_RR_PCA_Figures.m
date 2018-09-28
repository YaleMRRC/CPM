%This is a series of scripts that produce figures relevant to a PCA-based
%Ridge regression
%The following figures are:
% 1) Scree plot of variance explained by each PCA component
% 2) Boxplot of correlation coefficients across iterations
% 3) Actual versus predicted score scatter plots
% 4) Convert Edge Matrix to Connectivity matrix usable by BIS

%%
%%1) Scree plot of variance Explained by each component
clc;clear; close all
addpath /Users/danielbarron/GoogleDrive/NPC/CODE
addpath /Users/danielbarron/GoogleDrive/NPC/Transfer

list_phenotypes= {'wms','wais'};
title_phenotypes= {'WMS Variance Explained by Component','WAIS Variance Explained by Component'};

for x=1:length(list_phenotypes)
    load(sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/RR_6tasks_%s_PCA',list_phenotypes{x}));  
    explained = [];
    for j=1:length(results)
        for k=1:length(results(1).output_pca)
            explained(j,k,:) = results(j).output_pca(k).explained;
        end
    end
	explained = squeeze(mean(mean(explained,1),2));
    figure('units','normalized','outerposition',[0 0 1 1])
    pareto(explained)
    xlabel('Principal Component','FontSize',25)
    ylabel('Variance Explained (%)','FontSize',25)
    title(title_phenotypes{x},'FontSize',30); 
    print('-dpng',sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/RR_6tasks_%s_PCA_FIGURE_scree.png',list_phenotypes{x}))
    end
%%
%2) This produces a boxplot of the correlation coefficients across iterations.

clc;clear; close all
addpath /Users/danielbarron/GoogleDrive/NPC/CODE
addpath /Users/danielbarron/GoogleDrive/NPC/Transfer

list_phenotypes = {'wms','wais'};
%coefficients_both = zeros(length(list_phenotypes),10)
for x=1:length(list_phenotypes)
    load(sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/RR_6tasks_%s_PCA',list_phenotypes{x}));  
    for j=1:length(results)
    r_spearmann(j,:) = results(j).r_rank(:,1);
    r_pearson(j,:) = results(j).r_pearson(:,1);
    end
    %coefficients_both(:,x) = r_spearmann(:)
    coefficients_both(:,x) = r_pearson(:)
end

figure
boxplot([coefficients_both(:,1),coefficients_both(:,2)],'Notch','on','Labels',{'WMS Principle Component','WAIS Principle Component'})
ylim([0 0.7])
set(gca,'YTick',0.0:0.1:0.7)
ylabel('Pearson Correlation')
title('Pearson Correlation of Actual with Predicted Score')
print('-dpng',sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/Correlation_figures.png'))
   


%%
%3) Actual versus predicted score scatter plots, coloring individual dots with
%colors relevant to subject groups
clc;clear; close all
addpath /Users/danielbarron/GoogleDrive/NPC/CODE
addpath /Users/danielbarron/GoogleDrive/NPC/Transfer

list_phenotypes= {'wms','wais'};
title_phenotypes= {'WMS Actual versus Predicted Correlation','WAIS Actual versus Predicted Correlation'};

%Make the group assignments
%Trim down the file sizes to only those subjects included in the present
%analysis
delimiterIn = '\t';
groups = importdata('/Users/danielbarron/GoogleDrive/NPC/RAW/PhenotypeTXT/delete_group.txt',delimiterIn);
trimmed = groups(:,2);
delete = importdata('/Users/danielbarron/GoogleDrive/NPC/RAW/PhenotypeTXT/delete_leave236subs.txt',delimiterIn);
delete = nonzeros(delete(:,1));
trimmed(delete,:)=[];

%1=Control; 2=SCZ; 3=BPAD; 4=ADHD
[row1,col1] = find(trimmed==1);
[row2,col2] = find(trimmed==2); 
[row3,col3] = find(trimmed==3);
[row4,col4] = find(trimmed==4);

for x=1:length(list_phenotypes)
%Load the data
    load(sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/RR_6tasks_%s_PCA',list_phenotypes{x}));
    
       for i=1:size(results,2) %breaks by kfolds
        predicted_raw(i,:) = results(i).y;
        actual_raw(i,:) = results(i).pcascores_test_all(:,1); %%Potential error. Need to think thorugh the column selection here
       end
       actual = mean(actual_raw,1)';%%Potential error
       predicted = mean(predicted_raw,1)';
  
       hf=figure('units','normalized','outerposition',[0 0 1 1])
       hold on
        %Make the subplots 
        plot(actual(row1,col1),predicted(row1,col1),'k.','MarkerSize',35);%Normals are black
        plot(actual(row2,col2),predicted(row2,col2),'g.','MarkerSize',35);%Scz=green
        plot(actual(row3,col3),predicted(row3,col3),'b.','MarkerSize',35);%BPAD=blue
        plot(actual(row4,col4),predicted(row4,col4),'r.','MarkerSize',35);%ADHD=Red
        lsline;
        xlabel('Actual Principal Component Score','FontSize',25);
        ylabel('Predicted Principal Component Score','FontSize',25);
        title(title_phenotypes(x),'FontSize',30); 
        hold off
        print('-dpng',sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/RR_6tasks_%s_PCA_FIGURE_scatterplot.png',list_phenotypes{x}))
end
 %%   
%%4) Convert Edge Matrix to Connectivity matrix usable by BIS
clc;clear; close all
addpath /Users/danielbarron/GoogleDrive/NPC/CODE
addpath /Users/danielbarron/GoogleDrive/NPC/Transfer

list_phenotypes = {'wms','wais'};

for x=1:length(list_phenotypes)
%Load the data
    load(sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/RR_6tasks_%s_PCA',list_phenotypes{x}));  

    no_node = 268;
    ntasks = 6;
    no_kfolds = size(results,2);
    no_iterations = size(results(1).coef_total,2);

    aa = ones(no_node, no_node);
    aa_upp = triu(aa, 1);
    upp_id = find(aa_upp);

    MxMxN_matrix = zeros(no_node, no_node,ntasks,no_iterations,no_kfolds);
    for i=1:no_kfolds %breaks by kfolds
        for j=1:no_iterations %Breaks by iteration
            coef = results(i).coef_total(:,j); %this gives a 214668x1 vector
            for k=1:6 %breaks things down by task
                edge_vector = coef((1+(35778*(k-1))):(35778*k),:);
                %binarize
                edge_vector_bin = zeros(length(edge_vector),1);
                for l=1:length(edge_vector)
                    if edge_vector(l,:)>0
                        edge_vector_bin(l,:)=1;
                    end
                end   
                edge_vector_matrix = MxMxN_matrix(:,:,k);
                edge_vector_matrix(upp_id) = edge_vector_bin;
                edge_vector_matrix = edge_vector_matrix + edge_vector_matrix';
                MxMxN_matrix(:,:,k,j,i) =edge_vector_matrix;
            end
        end
    end
    final_nodes = mean(mean(mean(MxMxN_matrix,5),4),3);

    final_nodes_bin = zeros(268,268);
    for i=1:268
        for j=1:268
            if final_nodes(i,j)>=1.0
                final_nodes_bin(i,j) = 1; 
            end
        end
    end
    Percent_nodes = sum(sum(final_nodes_bin))/(268*268);
    fprintf('%0.2d of Possible Edges', Percent_nodes);
    dlmwrite(sprintf('/Users/danielbarron/GoogleDrive/NPC/Transfer/%s_FINAL_bin_matrix.txt',list_phenotypes{x}),final_nodes_bin, 'delimiter', '\t');
end
