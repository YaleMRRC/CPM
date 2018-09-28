%This loops through the ridgeCPM stuff
clc;clear; close all
addpath /data17/mri_group/daniel_data/npc/CODE
addpath /data17/mri_group/daniel_data/npc/Transfer

%Load data
load /data17/mri_group/daniel_data/npc/Transfer/all_mats_6tasks.mat

list_phenotypes = {'wms','wais'};
parpool(5)
for i=1:length(list_phenotypes)
    %Load phenotype
    phenotype = load(sprintf('/data17/mri_group/daniel_data/npc/Transfer/%s_raw',list_phenotypes{i}));
    all_behav = phenotype;
    parfor j=1:100
    fprintf('Run %.0d for %s',j,list_phenotypes{i}) ;
    seed = randi([1 10000]); thresh = 0.05;
   [results(j).q_s, results(j).r_pearson, results(j).r_rank, results(j).y, results(j).coef_total, results(j).coef0_total, results(j).lambda_total,results(j).output_pca,results(j).pcascores_test_all] = ridgeCPM_PCA(all_mats, all_behav,thresh,seed);
    end
    save(sprintf('/data17/mri_group/daniel_data/npc/Transfer/RR_6tasks_%s_PCA_100iterations',list_phenotypes{i}),'results');  
end


%% ExeFxn Tasks
%Define lists and variables to be used in CPM
list_matricies = {'bart','stopsignal','taskswitch'};
list_phenotypes= {'wais_lns','wais_mr','wais_voc','wms_ds','wms_ssp','wms_vr1ir','wms_vr2dr','hopkins_intsensitivity','hopkins_globalseverity','hopkins_somatization','hopkins_anxiety','hopkins_obscomp','hopkins_depression'};

all_mats = [];
all_behav = [];
for x=1:length(list_matricies)
    %Load and trim matrix file
    load(sprintf('/Users/danielbarron/GoogleDrive/NPC/MATRIX/%s_wnulls.mat',list_matricies{x}));
    M = squeeze(M);
    M(:,:,delete) = [];
    all_mats(:,:,:,x)=M;
end
    
for i=1:length(list_phenotypes)
    %Load and trim phenotype
    phenotype = load(sprintf('/Users/danielbarron/GoogleDrive/NPC/RAW/Phenotype_transformed/%s',list_phenotypes{i}));
%    phenotype(delete,:)=[];
    all_behav = phenotype;
    thresh = 0.05;
    [q_s, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total] = ridgeCPM(all_mats, all_behav, thresh);
    save(sprintf('/Users/danielbarron/GoogleDrive/NPC/CPM/RidgeRegression/RR_Exe_tasks_%s',list_phenotypes{i}),'q_s', 'r_pearson', 'r_rank', 'y', 'coef_total', 'coef0_total', 'lambda_total');
end


%%%MemFxn Tasks
%Define lists and variables to be used in CPM
list_matricies = {'pamenc','pamret','scap'};
list_phenotypes= {'wais_lns','wais_mr','wais_voc','wms_ds','wms_ssp','wms_vr1ir','wms_vr2dr','hopkins_intsensitivity','hopkins_globalseverity','hopkins_somatization','hopkins_anxiety','hopkins_obscomp','hopkins_depression'};

all_mats = [];
all_behav = [];
for x=1:length(list_matricies)
    %Load and trim matrix file
    load(sprintf('/Users/danielbarron/GoogleDrive/NPC/MATRIX/%s_wnulls.mat',list_matricies{x}));
    M = squeeze(M);
    M(:,:,delete) = [];
    all_mats(:,:,:,x)=M;
end
    
for i=1:length(list_phenotypes)
    %Load and trim phenotype
    phenotype = load(sprintf('/Users/danielbarron/GoogleDrive/NPC/RAW/Phenotype_transformed/%s',list_phenotypes{i}));
%    phenotype(delete,:)=[];
    all_behav = phenotype;
    thresh = 0.05;
    [q_s, r_pearson, r_rank, y, coef_total, coef0_total, lambda_total] = ridgeCPM(all_mats, all_behav, thresh);
    save(sprintf('/Users/danielbarron/GoogleDrive/NPC/CPM/RidgeRegression/RR_Mem_tasks_%s',list_phenotypes{i}),'q_s', 'r_pearson', 'r_rank', 'y', 'coef_total', 'coef0_total', 'lambda_total');
end
