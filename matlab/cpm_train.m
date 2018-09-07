function [r,p,pmask,mdl]=cpm_train(x,y,pthresh)
% call from predict_behavior
%posnegmask->pmask

% feature selection from all train subs
[r,p]=corr(x',y);
% pmask=p<pthresh & r>0;
% npmask=p<pthresh & r<0;

pmask=(+(r>0))-(+(r<0)); % BINARIZED
pmask=pmask.*(+(p<pthresh));

% get summary features for each train sub
for i=1:size(x,2)
%     summary_feature(i)=sum(x(pmask,i));
%    summary_feature(i)=mean(x(pmask,i))-mean(x(npmask,i));
    summary_feature(i)=mean(x(pmask>0,i))-mean(x(pmask<0,i));

end

% fit behavior to summary features
mdl=polyfit(summary_feature,y',1);
    
    