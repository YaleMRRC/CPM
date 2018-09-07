import numpy as np 
import scipy as sp
from scipy import io 
from scipy import signal as sg
#from matplotlib import pyplot as plt
#import seaborn as sns
import sklearn as skl
from sklearn import svm,cluster
import pandas as pd
import glob

def cosine_similarity(timeseries):
    """
    Function to calculate similarity between timeseries as a
    function of the angle of the complex representation
    Takes NxM matrix, where M = number of timeseries, and 
    N = number of timepoints
    Returns a matrix of size N x M x M
    """
    n_ts=timeseries.shape[1]
    n_tp=timeseries.shape[0]
    hilt = sg.hilbert(timeseries,axis=0)
    angles = np.angle(hilt)

    pw_diff=np.array([angles[v,:] - a for v in range(0,n_tp) for a in angles[v,:]])
    pw_diff=np.reshape(pw_diff,[n_tp,n_ts,n_ts])

    cos_sim=np.cos(pw_diff)

    return cos_sim


def calc_eigs(matrices,numevals="All"):
    """
    Takes NxMxM matrix and returns eigenvalues and eigenvectors
    """
    
    nvols,nrois,_=matrices.shape


    evals=np.zeros([nvols,nrois])
    evecs=np.zeros([nvols,nrois,nrois])
    evars=np.zeros([nvols,nrois])

    for volnum in range(0,nvols):

        eigs=sp.linalg.eigh(matrices[volnum,:,:])

        tevals=eigs[0]
        tevecs=eigs[1]

        evsort=np.argsort(tevals)
        tevals=tevals[evsort[::-1]]
        evals[volnum,:]=tevals
        evecs[volnum,:,:]=tevecs[:,evsort[-1::]]
        evars[volnum,:]=np.array([tevals[i]/np.sum(tevals,axis=None) for i in range(0,tevals.shape[0])])
 
        

        #evecs=np.array([evecs[i,:,evsort[i,:]] for i in range(0,len(evsort))])
        #evars=np.array([evals[i,:]/np.sum(evals[i,:],axis=None) for i in range(0,evals.shape[0])])


    opdict={}

    if numevals == 'All':
        opdict['EigVals']=evals
        opdict['EigVecs']=evecs
        opdict['EigVars']=evars

    else:
        opdict['EigVals']=evals[:,0:numevals]
        opdict['EigVecs']=evecs[:,:,0:numevals]
        opdict['EigVars']=evars[:,0:numevals]

    return opdict


def tsfilt(timeseries):
    """
    Demean and detrend input timeseries of one subject
    accepts array of size NxM, N = timepoints, M = timeseries
    """
    ts_detrend=sg.detrend(timeseries,axis=0)
    ts_demean=ts_detrend-ts_detrend.mean(axis=0)
    
    return ts_demean


def dotnorm(v1,v2):
    return  np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2)   

def indv_leida_mats(onesubdata,numeigs=1):
    """
    """

    filtered_data=tsfilt(onesubdata)

    cos_sim_data=cosine_similarity(filtered_data)

    opdict=calc_eigs(cos_sim_data,numevals=numeigs)

    evecs=opdict['EigVecs']
    tp,ts,numevec=evecs.shape
    fcd_list=[]
    for fcdi in range(0,numevec):
        evec=np.squeeze(evecs[:,:,fcdi])
        dns=np.array([dotnorm(e1, e2) for e1 in evec for e2 in evec])
        fcd_list.append(np.reshape(dns,[tp,tp]))


    opdict['FCD']=fcd_list

    return opdict

def group_leida_mats(subdict,nAs,numeigs):


    FCDdict={}
    for subid in subdict.keys():
        print('Processing Eigenvalues for sub: ',subid)
        ipdata=subdict[subid][:,0:nAs]
        eigdict=indv_leida_mats(ipdata,numeigs=numeigs)
        FCDdict[subid]=eigdict['FCD'][0]


    return FCDdict
        


def testcode():

    cossim=io.loadmat('../SW0065C_cossim.mat')
    cossimML=cossim['cossimmat']

    leML=io.loadmat('../SW0065C_leadeig.mat')
    leML=leML['Leading_Eig']

    fcd=io.loadmat('../SW0065C_FCD.mat')
    fcdML=fcd['FCD_eig']

    timeseries=io.loadmat('../Aging_data.mat')
    onesubdata=timeseries['SW0065C']
    onesubdata=onesubdata[:,:90]
    #onesubdata=np.vstack(onesubdata[5,:]).T

    filtered_data=tsfilt(onesubdata)

    cos_sim_data=cosine_similarity(filtered_data)

    opdict=calc_eigs(cos_sim_data,numevals=1)

    evecs=opdict['EigVecs']
    tp,ts,numevec=evecs.shape
    fcd_list=[]
    for fcdi in range(0,numevec):
        evec=np.squeeze(evecs[:,:,fcdi])
        dns=np.array([dotnorm(e1, e2) for e1 in evec for e2 in evec])
        fcd_list.append(np.reshape(dns,[tp,tp]))

    fcd=np.squeeze(np.array(fcd_list))
    #fcd[np.isclose(fcd,fcdML)]
    #fcd[np.isclose(fcd,fcdML)]=-fcd[np.isclose(fcd,fcdML)]
    plt.subplot(2,1,1)
    plt.imshow(fcd)
    plt.subplot(2,1,2)
    plt.imshow(fcdML)
    plt.show()


def kcluster_train(ipfeatures,nclust):
    est=cluster.KMeans(n_clusters=nclust)
    est.fit(ipfeatures)

    return est



def moody_data_split(rand=False):


    print("Loading data and creating subjects list")

    mats=glob.glob("/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST_LR/matrices/*REST*LR*_GSR*roimean.txt")

    pmatdf=pd.read_csv('/home/dmo39/pmatfilter.csv')
    substouse=list(map(str, pmatdf.Subject.values))

    matstouse=sorted([m for m in mats if any(s in m for s in substouse)])

    imn_file='/mnt/store2/mri_group/dave_data/code/LEiDA-master/LEiDA/netrank.mat'

    imnets_load=io.loadmat(imn_file)
    imnets=imnets_load['imnets']

    opname='/mnt/store2/mri_group/dave_data/code/LEiDA-master/LEiDA/NetParcelCorrMats_HCP.mat'


    if rand:
        print("You chose to randomize!")
        imnets=np.random.permutation(imnets)

        opimnets={}
        opimnets['imnets']=imnets

        io.savemat(imn_file.replace('.mat','_random.mat'),opimnets)



        opname=opname.replace('.mat','_randassign.mat')
 
    print("output file will be:", opname)

    split_ts(matstouse,imnets,opname)


def split_ts(ipmats, statelbls, savepath):
    #mats=glob.glob("/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST_LR/matrices/*REST*LR*_GSR*roimean.txt")

    #mats=glob.glob("/mnt/store1/*roimean.txt")

    matnames=[m.split('/')[-1] for m in ipmats]

    print("reading matrices")

    roidfs=[pd.read_csv(m,index_col=0,sep='\t').dropna(axis=1) for m in ipmats]

    data_dct={}

    print("starting splitting")
    for netnum in np.unique(statelbls):
        group_list=[]    

        for dfnum,roidf in enumerate(roidfs):

            print('Network ',netnum, 'Subject', dfnum)

            imnet=statelbls[:,dfnum]
            subdf=roidf[imnet==netnum]
            corrmat=subdf.corr().values
            group_list.append(corrmat)

        netmats=np.array(group_list)

        data_dct['Network'+str(netnum)]=netmats


    #io.savemat('/mnt/store2/mri_group/dave_data/code/LEiDA-master/LEiDA/NetParcelCorrMats_HCP.mat',data_dct)
    io.savemat(savepath,data_dct)

        
def svr_():

    imnets_load=io.loadmat('/mnt/store2/mri_group/dave_data/code/LEiDA-master/LEiDA/netrank.mat')
    imnets=imnets_load['imnets']

    pmatdf=pd.read_csv('/home/dmo39/pmatfilter.csv')
    excludevals=pmatdf[pmatdf.PMAT24_A_CR.isna()].index.values
    pmatdf_filt=pmatdf.dropna()
    pmatvals=pmatdf_filt.PMAT24_A_CR.values





    features=np.array([np.histogram(ims,bins=11)[0] for ims in imnets.T])

    svr=svm.SVR()

    param_grid = [
      {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
      {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
     ]

    clf = GridSearchCV(svr, param_grid, cv=5)

    clf.fit(features_filt_train,pmatvals_train)

    includevals=list(set(np.arange(0,865))-set(excludevals))

    features_filt=features[includevals,:]

    np.corrcoef(clf.predict(features_filt[401:,:]),pmatvals[401:])
