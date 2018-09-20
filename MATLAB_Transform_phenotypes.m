%Transform the phenotypes to be normal
%This file takes a list of phenotypes (which I've done for my dataset) and
%performs different transforms on it, producing a figure of a before/after
%histogram with skewness values. This does NOT automatically make data
%gaussian, but it's a first-show

clc;clear;close all
%Load the delete file
delimiterIn = '\t';
delete = importdata('/Users/danielbarron/GoogleDrive/NPC/RAW/PhenotypeTXT/delete_leave236subs.txt',delimiterIn);
delete = nonzeros(delete);

%Load the WAIS phenotypes
% list_pheno= {'wais_lns','wais_mr','wais_voc'};
% list_titles = {'Letter Number Sequencing','Matrix Reasoning','Vocab'};
% save_title = 'waisTransformed';

%Load the WMS phenotypes
list_pheno= {'wms_ds','wms_ssp','wms_vr1ir','wms_vr2dr',};
list_titles = {'Digit Span','Symbol Span','Visual Repro 1','Visual Repo 2'};
save_title = 'wmsTransformed'

%Load the Hopkins phenotypes
% list_pheno = {'hopkins_intsensitivity','hopkins_globalseverity','hopkins_somatization','hopkins_anxiety','hopkins_obscomp','hopkins_depression'};
% list_titles = {'Interpersonal Sensitivity','Global Severity','Somatization','Anxiety','Obsessive-Compulsive','Depression'};
% save_title = 'hopkinsTransformed'

figure('units','normalized','outerposition',[0 0 1 1])
 %WAIS is left-skewed, so use a squared transform
for i=1:length(list_pheno)
    phenotype = load(sprintf('/Users/danielbarron/GoogleDrive/NPC/RAW/Phenotype_matlab/%s.txt',list_pheno{i}));
    phenotype(delete) = [];
    %Plot raw data
    ha(i*2-1) = subplot(length(list_pheno),2,(i*2-1));
    hold on
    hist(phenotype(:))
    skew = skewness(phenotype);
    title(sprintf('RAW %s with Skewness=%.2f',list_titles{i},skew))
    hold off
   
    %Option 1: look at skewness and use either a square/cubed/etc if left
    %skewed and log if right-skewed
    if (skew < -.9)
        pheno = phenotype.^(4);
    elseif (skew < -0.7) && (skew > -0.9)
        pheno = phenotype.^(3);
    elseif (skew < -0.2) && (skew > -0.7)
        pheno = phenotype.^(2);
    elseif (skew > 1.0)
       pheno = log(phenotype/mean(phenotype)+0.2)+10; 
    elseif (skew > 0.2) && (skew < 1.0) 
       pheno = log(phenotype/mean(phenotype)+0.4)+10;
        %Log transform data and plot
    else
        pheno = phenotype;
    end
    %Option 2: Use BoxCox transform
    %pheno = boxcox(phenotype+1);
    
    %Option 3: do no transform
    %pheno = phenotype;

    skew_trans = skewness(pheno);
    ha(i*2) = subplot(length(list_pheno),2,i*2);
    hold on
    hist(pheno)
    title(sprintf('Transformed %s with Skewness=%.2f',list_titles{i},skew_trans))
    hold off
    dlmwrite(sprintf('/Users/danielbarron/GoogleDrive/NPC/RAW/Phenotype_transformed/%s_untransformed',list_pheno{i}),pheno, 'delimiter', '\t');
end

print('-dpng',sprintf('/Users/danielbarron/GoogleDrive/NPC/RAW/Phenotype_transformed/FIGURES/DIST_%s_logsq.png',save_title))

