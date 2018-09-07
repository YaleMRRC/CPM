import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from scipy import io
from matplotlib.colors import ListedColormap
from sklearn import linear_model

data_dir='C:/Users/david/Documents/Research/10simplerules_data_figs/data'
fig_dir='C:/Users/david/Documents/Research/10simplerules_data_figs/figs'


#### Funcs to use ####

def bias_plot():
    """
    Bias plot for 10 simple rules figure
    """

    biasstuff=io.loadmat(data_dir+'/cv_compare_tot_norm.mat',mat_dtype=True)
    bias_data=biasstuff['rsn_tot'][0]
    #cvtypes=['k2','k5','k10','loo','external']
    #cvtypes_plots=['Split Half', '5 Fold', '10 Fold', 'LOO', 'External']

    cvtypes=['k2','k5','k10','loo']
    cvtypes_plots=['Split-half', '5-fold', '10-fold', 'LOO']


    Rpos=np.concatenate([bias_data[0][cvt][:,0] for cvt in cvtypes])
    Rmse=np.concatenate([bias_data[0][cvt][:,4] for cvt in cvtypes])
    labels=np.concatenate([np.repeat(cv,200) for cv in cvtypes_plots])

    dfarr=np.concatenate([np.vstack(labels).T,np.vstack(Rpos).T,np.vstack(Rmse).T]).T

    #n1="R Squared (Pearsons)"
    #n2="R Squared (MSE)"
    n1="MSE(observed,yhat)"
    n2="MSE(observed,pred)"

    ipdata=pd.DataFrame(dfarr,columns=['Labels',n1,n2])

    ipdata[n1]=ipdata[n1].astype('float')
    ipdata[n2]=ipdata[n2].astype('float')

    ipdata[n1]=ipdata[n1]**2
    ipdata[n2]=ipdata[n2]

    ipdata['X']=np.concatenate([np.linspace(0,0.3,200) for n in range(0,len(cvtypes_plots))])
    ipdata['Y']=np.concatenate([np.linspace(0,0.3,200) for n in range(0,len(cvtypes_plots))])


    g=sns.FacetGrid(ipdata,col="Labels",col_wrap=2,hue="Labels",palette='vlag',despine=True)
    axes = g.axes.flatten()
    axes[0].set_title(r"$Split-half$")
    axes[1].set_title(r"$5-fold$")
    axes[2].set_title(r"$10-fold$")
    axes[3].set_title(r"$LOO$")
    g=g.map(plt.plot,"X","Y")
    g=g.map(plt.scatter,n1,n2).set_axis_labels(r'$MSE(observed,\hat y)$', r'$MSE(observed,pred)$')


    axes[0].set_title("Split-half")
    axes[1].set_title("5-fold")
    axes[2].set_title("10-fold")
    axes[3].set_title("LOO")

    #plt.show()

    print(fig_dir+'/biasfig.tiff')
    plt.savefig(fig_dir+'/biasfig.tiff', dpi=300, facecolor='w', edgecolor='w',
         orientation='portrait', papertype=None, format=None,
         transparent=False, bbox_inches=None, pad_inches=0.1,
         frameon=None)   

#    sns.lmplot(x="PearsonsR", y="mseR",col="Labels", hue="Labels", data=ipdata,col_wrap=3,ci=None,fit_reg=True)
#    plt.show()

    


    # for i,cv in enumerate(ipdata.Labels.unique()):
    #     g=plt.subplot(2,3,i+1)
    #     print(cv)
    #     subdf=ipdata[ipdata.Labels == cv]
    #     regr = linear_model.LinearRegression()
    #     regr.fit(np.vstack(subdf.PearsonsR.values), subdf.mseR.values)
    #     #print(' '.join([cv,str(regr.coef_[0]),str(regr.intercept_)]))
    #     g=sns.FacetGrid(subdf)
    #     g=g.map(plt.scatter,"PearsonsR","mseR")
    #     plt.plot([0,0.55],[0,0.55],color='k')
    # plt.show()

    

def nsubs_plot():
    # Ntrainsubs plot
    nsubsdata=io.loadmat(data_dir+'/norm_nsubtrain.mat')
    nsubsdata=nsubsdata['behav_struct_nsubs']
    
    train_lbls=['train'+str(num) for num in np.arange(25,401,25)]    
    col_lbls=[str(num) for num in np.arange(25,401,25)]

    #rpos=np.array([nsubsdata['res_struct_nsubs'][0][tl][0][:,1] for tl in train_lbls]).T
    
    # Aggregate predicted behavior-
    pred_data_dict={cvt:nsubsdata[0][cvt][0]['predbehavpos'][0][0] for cvt in train_lbls}
    real_data_dict={cvt:nsubsdata[0]['testbehav'][0] for cvt in train_lbls}
    

    cllct=[]
    for cvt in train_lbls:
    
        mse=np.mean(np.array(pred_data_dict[cvt]-real_data_dict[cvt])**2,axis=1)
        cllct.append(mse)

    mse_tot=np.array(cllct)
    
    ipdata=pd.DataFrame(mse_tot.T,columns=col_lbls)


    #rposdf=pd.DataFrame(rpos,columns=col_lbls)
    
    train_subs_plot(ipdata, fig_dir+'/nsubs_mse.tiff')


def sametrain_run():
    # Const train size plot
    ctrain_data=io.loadmat(data_dir+'/sametrain_norm.mat')
    #cvdata=ctrain_data['res_struct_norm_st']
    cvdata=ctrain_data['behav_struct_norm_st']
    #cvtypes=['k2','k5','k10','loo','external']
    #cvtypes_plots=['Split Half', '5 Fold', '10 Fold', 'LOO', 'External']
 
    cvtypes=['k2','k5','k10','loo']
    cvtypes_plots=['Split-half', '5-fold', '10-fold', 'LOO']

    #rpos=np.array([ctrain_data[0][cvt][0][:,0] for cvt in cvtypes]).T

    #rposdf=pd.DataFrame(rpos,columns=cvtypes_plots)

    pncpmat=io.loadmat(data_dir+"/pncpmats.mat",variable_names=['pmats_pnc'])
    pncpmat=pncpmat['pmats_pnc'][0]
    pncpmat=(pncpmat-pncpmat.mean())/pncpmat.std()

    # Aggregate predicted behavior-
    pred_data_dict={cvt:cvdata[0][cvt][0]['predbehavpos'][0][0] for cvt in cvtypes}
    real_data_dict={cvt:cvdata[0][cvt][0]['testbehav'][0][0] for cvt in cvtypes}
    real_data_dict['external']=np.repeat(np.vstack(pncpmat).T,100,axis=0)


    cllct=[]
    for cvt in cvtypes:
    
        mse=np.mean(np.array(pred_data_dict[cvt]-real_data_dict[cvt])**2,axis=1)
        cllct.append(mse)

    mse_tot=np.array(cllct)
    
    ipdata=pd.DataFrame(mse_tot.T,columns=cvtypes_plots)

    ipdata=ipdata[['Split-half', '5-fold', '10-fold', 'LOO']]

    const_trainsize_plot(ipdata,fig_dir+'/sametrain_mse.tiff')

def mseplot():
    # Load Data
    cvdata=io.loadmat(data_dir+'/norm_cvcomparison_cuda.mat')
    cvdata=cvdata['behav_struct_norm']
    cvdata_extra=io.loadmat(data_dir+'/norm_cvcomparison_shred.mat')
    cvdata_extra=cvdata_extra['behav_struct_norm']
    pncpmat=io.loadmat(data_dir+"/pncpmats.mat",variable_names=['pmats_pnc'])
    pncpmat=pncpmat['pmats_pnc'][0]
    pncpmat=(pncpmat-pncpmat.mean())/pncpmat.std()
    # Set labels
    cvtypes=['k2','k5','k10','loo','external']
    cvtypes_plots=['Split-half', '5-fold', '10-fold', 'LOO', 'External']
    cvtypes_int=['k2','k5','k10','loo']
    
    # Aggregate predicted behavior-
    pred_data_dict={cvt:cvdata[0][cvt][0]['predbehavpos'][0][0] for cvt in cvtypes}
    real_data_dict={cvt:cvdata[0][cvt][0]['testbehav'][0][0] for cvt in cvtypes_int}
    real_data_dict['external']=np.repeat(np.vstack(pncpmat).T,100,axis=0)


    cllct=[]
    for cvt in cvtypes:

        mse=np.mean(np.array(pred_data_dict[cvt]-real_data_dict[cvt])**2,axis=1)
        cllct.append(mse)

    mse_tot=np.array(cllct)
    
    ipdata=pd.DataFrame(mse_tot.T,columns=cvtypes_plots)

    ipdata=ipdata[['Split-half', '5-fold', '10-fold', 'LOO']]

    cv_effect_plot(ipdata,fig_dir+'/normcvMSE.tiff')
    #cv_effect_plot_sepscales(ipdata,'../normcvMSE_splitscale.jpg')

def cv_effect_plot(ipdata,opname):
    plt.clf()
    sns.set(style='ticks')
    #sns.boxplot(data=ipdata[['Split Half', '5 Fold', '10 Fold', 'LOO']],palette='vlag')
    sns.boxplot(data=ipdata,palette='vlag')
    
    #sns.violinplot(data=rposdf,palette='vlag',inner='quartile')
    sns.despine(offset=10, trim=True)
    plt.title('Variable train size',fontsize='large')
    plt.ylabel('MSE',weight='bold')
    plt.xlabel('Cross-validation method',weight='bold')
    plt.tight_layout()
    #plt.savefig('./cv_pred_200.jpg', dpi=None, facecolor='w', edgecolor='w',
    #    orientation='portrait', papertype=None, format=None,
    #    transparent=False, bbox_inches=None, pad_inches=0.1,
    #    frameon=None)
    plt.savefig(opname, dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)



def train_subs_plot(ipdata,opname):
    plt.clf()
    sns.set(style='ticks')
    sns.boxplot(data=ipdata,palette='vlag')
    #sns.violinplot(data=rposdf,palette='vlag',inner='quartile')
    sns.despine(offset=10, trim=True)
    plt.title('Effect of number of training individuals on prediction',fontsize='large')
    plt.ylabel('MSE',weight='bold')
    plt.xlabel('Number of training individuals',weight='bold')
    plt.tight_layout()
    plt.savefig(opname, dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)


def const_trainsize_plot(ipdata,opname):
    plt.clf()
    sns.set(style='ticks')
    sns.boxplot(data=ipdata,palette='vlag')
    #sns.violinplot(data=rposdf,palette='vlag',inner='quartile')
    sns.despine(offset=10, trim=True)
    plt.title('Constant train size (n=180)',fontsize='large')
    plt.ylabel('MSE',weight='bold')
    plt.xlabel('Cross-validation method',weight='bold')
    plt.tight_layout()
    plt.savefig(opname, dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)


#### WIP ####

def run_cv_effect(fpath,structname,opname):
    cvdata=io.loadmat(fpath)

    cvdata=cvdata[structname]

    cvtypes=['k2','k5','k10','loo','external']
    
    cvtypes_plots=['Split Half', '5 Fold', '10 Fold', 'LOO', 'External']
    
    rpos=np.array([cvdata[0][cvt][0][:,0] for cvt in cvtypes]).T

    ipdata=pd.DataFrame(rpos,columns=cvtypes_plots)
    #ipdata=ipdata[['Split Half', '5 Fold', '10 Fold', 'LOO']]

    cv_effect_plot(ipdata,opname)



def prior_dataproc():
    # Initial CV comparison plot
    cvdata=io.loadmat('./iteres.mat')
    cvdata=cvdata['res_struct']
    #cvdata_extra=io.loadmat('./iteres_extra.mat')
    #cvdata_extra=cvdata_extra['res_struct_extra']
    cvtypes=['k2','k5','k10','loo','external']
    cvtypes_plots=['Split Half', '5 Fold', '10 Fold', 'LOO', 'External']
    rpos=np.array([cvdata[0][cvt][0][:,0] for cvt in cvtypes]).T
    #rpos_extra=np.array([cvdata_extra[0][cvt][0][:,0] for cvt in cvtypes]).T
    #rposdf=pd.DataFrame(np.concatenate([rpos,rpos_extra]),columns=cvtypes_plots)
    rposdf=pd.DataFrame(np.concatenate([rpos,rpos_extra]),columns=cvtypes_plots)
    cv_effect_plot()





def quickcorrplot():
    # Initial CV comparison plot
    #cvdata=io.loadmat('../norm_cvcomparison_cuda.mat')
    #cvdata=cvdata['res_struct_norm']
    #cvdata_extra=io.loadmat('../norm_cvcomparison_shred.mat')
    #cvdata_extra=cvdata_extra['res_struct_norm']
    
    cvdata=io.loadmat(data_dir+'/LEresstruct.mat')
    cvdata=cvdata['res_struct']

    #cvtypes=['k2','k5','k10','loo','external']
    cvtypes=['k2','k5','k10','loo']
    
    #cvtypes_plots=['Split Half', '5 Fold', '10 Fold', 'LOO', 'External']
    cvtypes_plots=['Split-half', '5-fold', '10-fold', 'LOO']
    
    rpos=np.array([cvdata[0][cvt][0][:,0] for cvt in cvtypes]).T
    #rpos_extra=np.array([cvdata_extra[0][cvt][0][:,0] for cvt in cvtypes]).T
    #rposdf=pd.DataFrame(np.concatenate([rpos,rpos_extra]),columns=cvtypes_plots)
    rposdf=pd.DataFrame(rpos,columns=cvtypes_plots)

    cv_effect_plot(rposdf,fig_dir+'/LEcvR.jpg')


def cv_effect_plot_sepscales(ipdata,opname):
    plt.clf()


    #fig.set_ylabel('R',weight='bold')
    #fig.set_xlabel('Cross Validation Method',weight='bold')
    #my_cmap=ListedColormap(sns.color_palette('vlag').as_hex())
    my_cmap=sns.color_palette('vlag').as_hex()[0:5]

    f,a1=plt.subplots()
    sns.set_palette("vlag")

    props = dict(widths=0.7,patch_artist=True, medianprops=dict(color="black"))

    b1=a1.boxplot(ipdata[['Split Half', '5 Fold', '10 Fold', 'LOO']].values,positions=[0,1,2,3],**props)

    a2=a1.twinx()
    b2=a2.boxplot(ipdata[['External']].values,positions=[4],**props)
    
    a1.set_xlim(-0.5,4.5)
    a1.set_xticks(range(len(ipdata.columns)))
    a1.set_xticklabels(ipdata.columns)

    for patch, color in zip(b1['boxes']+b2['boxes'], my_cmap):
        patch.set_facecolor(color)

    sns.despine(offset=10, trim=True,right=False)
    sns.set(style='ticks')

    plt.suptitle('Variable train size',fontsize='large')
    a1.set_ylabel('MSE',weight='bold')
    a1.set_xlabel('Cross-validation method',weight='bold')

    plt.tight_layout()

    plt.savefig(opname, dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)    
