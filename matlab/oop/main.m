%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  N is the number of subjects
%  T is the number of tasks
%  x: 268*268*N*T
%  y: N*1
%  d: number of elements in lower triangular part of 268*268
%  group: 1*N of type subject
%  subject: a class of 9 task-based connectome and single label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = main()
%     data = csvread('../data/gender.csv');
    dataset = "hcp.50"; % LDA on UCLA + ages on 3 bins for HCP
    x = load('../data.50/HCP900_rest_n50.mat');
    y = load('../data.50/HCP900_PMAT24_A_CR_n50.mat');
    x= x.HCP900_rest_n50; % makes sure x is 268*268*N
    y=y.HCP900_PMAT24_A_CR_n50; % makes sure y is N*1
    N = size(x,3);
    g = buildGroup(x,dataset,'none'); % mask=false, Bins
    options = [];
    options.thresh=0.05;
    options.seed = randi([1 10000]);
    options.k = 2;
    options.phenotype = phenotype('behav',y);
    options.diagnosis = zeros(N,1);
    m = rcpm(g,options);
    m.run();
    m.evaluate();
end

function g = buildGroup(x,dataset,mask)
N =size(x,3);
subjects(1,N) = subject(N);
for i=1:N
    subjects(i) = subject(x(:,:,i,:),i,dataset,mask);
end
g = group(subjects);
end
