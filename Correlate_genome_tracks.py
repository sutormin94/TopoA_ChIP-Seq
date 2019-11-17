###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to calculate and plot correlation matrix of a set of genome tracks (wig).
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import scipy
import scipy.cluster.hierarchy as sch

#Dictionary of replicas 
#'Track name' : 'Path to wig file'
Dict_of_tracks_1={'ParC_1' : "F:\Sayyed_TopoIV_data\Cov_depth\ParC_1_FE.wig",
                  'ParE_1' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_1_FE.wig",
                  'ParE_2' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_2_FE.wig",
                  'ParE_G1' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_G1_FE.wig",
                  'ParE_S20min' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_S20min_FE.wig",
                  'ParE_S40min' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_S40min_FE.wig",  
                  'ParE_G2' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_G2_FE.wig", 
                  'NorflIP_ParC' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParC_1_FE.wig", 
                  'NorflIP_ParE' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_1_FE.wig",
                  'NorflIP_ParE' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_2_FE.wig",
                  'Cfx_10mkM_2' : "F:\TopoIV_Topo-Seq\Fold_enrichment\Cfx_10mkM_2_FE.wig",
                  'S83L_Cfx_10mkM' : "F:\TopoIV_Topo-Seq\Fold_enrichment\S83L_Cfx_10mkM_FE.wig",
                  'S83L_Cfx_100mkM' : "F:\TopoIV_Topo-Seq\Fold_enrichment\S83L_Cfx_100mkM_FE.wig"}
Dict_of_tracks={'TopA_CTD-Rif-_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_1_FE.wig",
                'TopA_CTD-Rif-_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_2_FE.wig",
                'TopA_CTD-Rif-_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_3_FE.wig",
                'TopA_CTD-Rif+_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_plus_1_FE.wig",
                'TopA_CTD-Rif+_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_plus_2_FE.wig",
                'TopA_CTD-Rif+_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_plus_3_FE.wig",
                'TopA_CTD+Rif-_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_minus_1_FE.wig",
                'TopA_CTD+Rif-_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_minus_2_FE.wig",
                'TopA_CTD+Rif-_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_minus_3_FE.wig",
                'TopA_CTD+Rif+_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_plus_2_FE.wig",
                'TopA_CTD+Rif+_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_plus_3_FE.wig",}
#Output path for the corralation matrix.
Outpath="C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_all_conditions_replicas_FE_correlation_matrix.png"
Outpath_clusterized="C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_all_conditions_replicas_FE_correlation_matrix_clusterized.png"


#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    return NE_values

#Contains data of all replicas in separate arrays.
dict_of_replicas={}
samples_names_array=[]
for replica_name, replica_path in Dict_of_tracks.items():
    samples_names_array.append(replica_name)
    dict_of_replicas[replica_name]=wig_parsing(replica_path)

Replicas_dataframe=pd.DataFrame(dict_of_replicas, columns=samples_names_array)
print(Replicas_dataframe.head(10))


#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath):
    fig=plt.figure(figsize=(8,8), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    cax=ax1.imshow(df.corr(method=cor_method), interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00], shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()
    plt.close()
    return

#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#######

def Clustering(df):
    X = df.corr().values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
    df = df.reindex_axis(columns, axis=1)
    return df

correlation_matrix(Replicas_dataframe, 'pearson', 'Correlation of samples', Outpath)
Replicas_dataframe_clusterized=Clustering(Replicas_dataframe)
correlation_matrix(Replicas_dataframe_clusterized, 'pearson', 'Correlation of samples clusterized', Outpath_clusterized)