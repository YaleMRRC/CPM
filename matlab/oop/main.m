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
    dataset = "hcp.175"; % LDA on UCLA + ages on 3 bins for HCP
    x = rand(268,268,175,6);%load('all_mats.mat');
    y = randi(30,175,1);%load('IQ.mat');
    g = buildGroup(x,dataset,'none'); % mask=false, Bins
    options = [];
    options.thresh=0.05;
    options.seed = randi([1 10000]);
    options.k = 2;
    options.phenotype = phenotype('behav',y);
    options.diagnosis = randi(2,175,1);
    m = cpm(g,options);
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
