import numpy as np 
import scipy as sp
import pandas as pd
#from matplotlib import pyplot as plt
#import seaborn as sns
import glob
from scipy import stats
import random


def read_mats(iplist):

    x=[pd.read_csv(m,sep='\t',header=None) for m in iplist]
    x=[df.dropna(axis=1).values for df in x]
    ipmats=np.stack(x,axis=2)

    return ipmats



def train_cpm(ipmat,pheno):

    """
    Accepts input matrices and pheno data
    Returns model
    """

    cc=[stats.pearsonr(pheno,im) for im in ipmat]



    rmat=np.array([c[0] for c in cc])
    pmat=np.array([c[1] for c in cc])
    rmat=np.reshape(rmat,[268,268])
    pmat=np.reshape(pmat,[268,268])
    posedges=(rmat > 0) & (pmat < 0.01)
    posedges=posedges.astype(int)
    negedges=(rmat < 0) & (pmat < 0.01)
    negedges=negedges.astype(int)
    pe=ipmat[posedges.flatten().astype(bool),:]
    ne=ipmat[negedges.flatten().astype(bool),:]
    pe=pe.sum(axis=0)/2
    ne=ne.sum(axis=0)/2


    if np.sum(pe) != 0:
        fit_pos=np.polyfit(pe,pheno,1)
    else:
        fit_pos=[]

    if np.sum(ne) != 0:
        fit_neg=np.polyfit(ne,pheno,1)
    else:
        fit_neg=[]

    return fit_pos,fit_neg,posedges,negedges


def pairwise_corr(X,Y):

    #A=A.astype('float64')

    n=X.shape[0]
    
    xbar=X.mean(1)
    ybar=Y.mean(1)

    xminusxbar=X-np.vstack(xbar)
    yminusybar=Y-np.vstack(ybar)

    xminusxbarsquare=xminusxbar**2
    yminusybarsquare=yminusybar**2

    numer=np.dot(xminusxbar.T,yminusybar)
    denom=np.vstack(np.sqrt(xminusxbarsquare.sum(1))*np.sqrt(yminusybarsquare.sum(1)))

    return numer/denom


def run_validate(ipmats,pheno,cvtype):

    numsubs=ipmats.shape[2]
    ipmats=np.reshape(ipmats,[-1,numsubs])

    

    if cvtype == 'LOO':
        behav_pred_pos=np.zeros([numsubs])
        behav_pred_neg=np.zeros([numsubs])
        for loo in range(0,numsubs):

            print("Running LOO, sub no:",loo)
      
            train_mats=np.delete(ipmats,[loo],axis=1)
            train_pheno=np.delete(pheno,[loo],axis=0)
            
            test_mat=ipmats[:,loo]
            test_phenp=pheno[loo]

            pos_fit,neg_fit,posedges,negedges=train_cpm(train_mats,train_pheno)

            pe=np.sum(test_mat[posedges.flatten().astype(bool)])/2
            ne=np.sum(test_mat[negedges.flatten().astype(bool)])/2

            if len(pos_fit) > 0:
                behav_pred_pos[loo]=pos_fit[0]*pe + pos_fit[1]
            else:
                behav_pred_pos[loo]='nan'

            if len(neg_fit) > 0:
               behav_pred_neg[loo]=neg_fit[0]*ne + neg_fit[1]
            else:
                behav_pred_neg[loo]='nan'

        
        Rpos=stats.pearsonr(behav_pred_pos,pheno)[0]
        Rneg=stats.pearsonr(behav_pred_neg,pheno)[0]

        return Rpos,Rneg


    elif cvtype == '5k':
        bp,bn,ba=kfold_cpm(ipmats,pheno,numsubs,5)



        ccp=np.array([stats.pearsonr(bp[i,:],ba[i,:]) for i in range(0,5)])
        Rpos_mean=ccp.mean(axis=0)[0]

        ccn=np.array([stats.pearsonr(bn[i,:],ba[i,:]) for i in range(0,5)])
        Rneg_mean=ccn.mean(axis=0)[0]



    elif cvtype == '10k':
        bp,bn,ba=kfold_cpm(ipmats,pheno,numsubs,10)


        ccp=np.array([stats.pearsonr(bp[i,:],ba[i,:]) for i in range(0,10)])
        Rpos_mean=ccp.mean(axis=0)[0]

        ccn=np.array([stats.pearsonr(bn[i,:],ba[i,:]) for i in range(0,10)])
        Rneg_mean=ccn.mean(axis=0)[0]



    elif cvtype == 'splithalf':
        bp,bn,ba=kfold_cpm(ipmats,pheno,numsubs,2)

        ccp=np.array([stats.pearsonr(bp[i,:],ba[i,:]) for i in range(0,2)])
        Rpos_mean=ccp.mean(axis=0)[0]

        ccn=np.array([stats.pearsonr(bn[i,:],ba[i,:]) for i in range(0,2)])
        Rneg_mean=ccn.mean(axis=0)[0]


    else:
        raise Exception('cvtype must be LOO, 5k, 10k, or splithalf')


    return Rpos_mean,Rneg_mean
    


def kfold_cpm(ipmats,pheno,numsubs,k):
    randinds=np.arange(0,numsubs)
    random.shuffle(randinds)

    samplesize=int(np.floor(float(numsubs)/k))

    behav_pred_pos=np.zeros([k,samplesize])
    behav_pred_neg=np.zeros([k,samplesize])

    behav_actual=np.zeros([k,samplesize])

    for fold in range(0,k):
        print("Running fold:",fold+1)
        si=fold*samplesize
        fi=(fold+1)*samplesize


        if fold != k-1:
            testinds=randinds[si:fi]
        else:
            testinds=randinds[si:]

        traininds=~np.isin(randinds,testinds)
        
        trainmats=ipmats[:,traininds]
        trainpheno=pheno[traininds]
 
        testmats=ipmats[:,testinds]
        testpheno=pheno[testinds]

        behav_actual[fold,:]=testpheno


        pos_fit,neg_fit,posedges,negedges=train_cpm(trainmats,trainpheno)

        pe=np.sum(testmats[posedges.flatten().astype(bool),:], axis=0)/2
        ne=np.sum(testmats[negedges.flatten().astype(bool),:], axis=0)/2


        if len(pos_fit) > 0:
            behav_pred_pos[fold,:]=pos_fit[0]*pe + pos_fit[1]
        else:
            behav_pred_pos[fold,:]='nan'

        if len(neg_fit) > 0:
            behav_pred_neg[fold,:]=neg_fit[0]*ne + neg_fit[1]
        else:
            behav_pred_neg[fold,:]='nan'

    return behav_pred_pos,behav_pred_neg,behav_actual


def sample_500(ipmats,pheno,cvtype):

    numsubs=ipmats.shape[2]

    randinds=np.arange(0,numsubs)
    random.shuffle(randinds)

    randinds500=randinds[:500]

    ipmats_rand=ipmats[:,:,randinds500]
    pheno_rand=pheno[randinds500]

    opdict={}

    Rpos_loo,Rneg_loo=run_validate(ipmats_rand,pheno_rand,'LOO')
    
    Rpos_2k,Rneg_2k=run_validate(ipmats_rand,pheno_rand,'splithalf')

    Rpos_5k,Rneg_5k=run_validate(ipmats_rand,pheno_rand,'5k')

    Rpos_10k,Rneg_10k=run_validate(ipmats_rand,pheno_rand,'10k')

    opdict['LOO_Rpos'] = Rpos_loo
    opdict['LOO_Rneg'] = Rneg_loo
    opdict['2k_Rpos'] = Rpos_2k
    opdict['2k_Rneg'] = Rneg_2k
    opdict['5k_Rpos'] = Rpos_5k
    opdict['5k_Rneg'] = Rneg_5k
    opdict['10k_Rpos'] = Rpos_10k
    opdict['10k_Rneg'] = Rneg_10k
    opdict['Sample_Indices']=randinds500

    return opdict



def shred_data_run_hcp():
    mats=glob.glob('*WM*LR*_GSR*.txt')
    mats=list(sorted(mats))
    pheno=pd.read_csv('unrestricted_dustin_6_21_2018_20_47_17.csv')
    subs=[m.split('_')[0] for m in mats]
    pfilesubs=set(pheno.Subject)
    subs=[int(s) for s in subs]
    subset=set(subs)
    usesubs=list(pfilesubs.intersection(subset))
    usesubs=list(map(str,usesubs))
    usesubs=sorted(usesubs)
    iplist=[m for m in mats if any([u in m for u in usesubs])]
    x=[pd.read_csv(m,sep='\t',header=None) for m in iplist]
    x=[df.dropna(axis=1).values for df in x]
    ipmats=np.stack(x,axis=2)
    phenofilt=pheno[pheno.Subject.isin(usesubs)]
    pmatvals=phenofilt['PMAT24_A_CR'].values

    return ipmats,pmatvals,usesubs


def shred_data_run_pnc():
    iplist=sorted(glob.glob('*matrix.txt'))
    pheno=pd.read_csv('phenotypes_info.csv')
    df_filter=pheno[['SUBJID','pmat_cr']]
    df_filter=df_filter.dropna()
    subs_mats=[i.split('_')[0] for i in iplist]
    subs_pheno=list(map(str,df_filter.SUBJID.unique()))
    subs_pheno=[sp.split('.')[0] for sp in subs_pheno]
    substouse=sorted(list(set(subs_mats) & set(subs_pheno)))

    iplist=[ip for ip in iplist if any([s in ip for s in substouse])]

    mats=[pd.read_csv(m,sep='\t',header=None).dropna(axis=1).values for m in iplist]
    ipmats=np.stack(mats,axis=2)
    pmatvals=df_filter.sort_values('SUBJID').pmat_cr.values
    

    return ipmats,pmatvals,substouse