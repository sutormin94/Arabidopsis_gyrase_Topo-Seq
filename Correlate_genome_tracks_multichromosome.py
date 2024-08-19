###############################################
##Dmitry Sutormin, 2024##
##ChIP-Seq analysis##

####
# Script calculates and plots a correlation matrix of a set of genome tracks (WIG).
# Handles multichromosome datasets.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import scipy.stats as sci
import scipy.cluster.hierarchy as sch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#Path to folder with wig files.
PWD_wig="/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/NUMT_masked/WIG/"

#Input: Continuous data (e.g. EcTopoI fold enrichment) (WIG).
#Dictionary of replicas 
#'Track name' : 'Path to wig file'
Dict_of_tracks={'TAIR10.1_1' :        PWD_wig + 'TAIR10.1_ChIP_Seq_NUMT_masked_FE_1.wig',
                'TAIR10.1_2' :        PWD_wig + 'TAIR10.1_ChIP_Seq_NUMT_masked_FE_2.wig',
                }

# Dictionoary with chromosomal regions to be masked: {chr_name : [[start, stop], ], }
Dict_of_deletions={}

# Correlation method: 'pearson' or 'spearman' for Pearson and Spearman correlation coefficients, respectively.
Correlation_method='pearson'

#Path to folder with output files.
PWD_out="/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/NUMT_masked/WIG/"


#######
##Parses WIG file.
#######

def wig_parsing(wigfile):
    
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    Dict_of_chromosomes_data={}
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0]=='fixedStep':
            chrom_name=line[1].split('=')[1]
            Dict_of_chromosomes_data[chrom_name]=[]
        if line[0] not in ['track', 'fixedStep']:
            Dict_of_chromosomes_data[chrom_name].append(float(line[0]))
    wigin.close()
    
    return Dict_of_chromosomes_data

#######
##Read, combine data into dataframe.
#######

def read_combine(Dict_of_tracks):
    
    #Contains data of all replicas in separate arrays.
    dict_of_replicas={}
    samples_names_array=[]
    i=0
    for replica_name, replica_path in Dict_of_tracks.items():
        i+=1
        print('Progress: ' + str(i) + '/' + str(len(Dict_of_tracks)))
        samples_names_array.append(replica_name)
        dict_of_replicas[replica_name]=wig_parsing(replica_path)
    
    return dict_of_replicas

#########
##Mask regions (deletions and multiplicated genes).
#########

def mask_regions(Rep_chrom_dict, dict_of_deletions):
    
    if len(dict_of_deletions)>0:
        
        Rep_chrom_dict_masked={}
        
        for rep_name, chrom_dict in Rep_chrom_dict.items():
            
            chrom_dict_masked={}
            
            for chrom_name, chrom_data in chrom_dict.items():
                
                if chrom_name in dict_of_deletions:
                    
                    #Remove items falling into delelted or masked regions.
                    regs_to_mask=dict_of_deletions[chrom_name]
                    maska=[0]*len(chrom_data)
                    
                    if len(regs_to_mask)>0:
                        for reg_to_mask in regs_to_mask:
                            reg_start=reg_to_mask[0]
                            reg_end=reg_to_mask[0]
                            maska[reg_start:reg_end]=[1]*(reg_end-reg_start)
                    
                    chrom_data_masked=[]
                    
                    for i in range(len(chrom_data)): 
                        if maska[i]==0:
                            chrom_data_masked.append(chrom_data[i])
                            
                    chrom_dict_masked[chrom_name]=chrom_data_masked
                            
                else:
                    
                    chrom_dict_masked[chrom_name]=chrom_data
                      
                print(f'Length of {chrom_name} from {rep_name} sample before masking: {len(chrom_dict[chrom_name])}')
                print(f'Length of {chrom_name} from {rep_name} sample after masking: {len(chrom_dict_masked[chrom_name])}')
        
            Rep_chrom_dict_masked[rep_name]=chrom_dict_masked
            
        return Rep_chrom_dict_masked
    
    else:
        print('No regions to mask.')
    
        return Rep_chrom_dict

#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def make_correlation_matrix_plot(Rep_chrom_dict_masked, cor_method, title, outpath_folder, file_name):
    
    # Create empty correlation matrix.
    rep_chrom_names=[]
    
    for rep_name, chrom_dict in Rep_chrom_dict_masked.items():
        
        chrom_name_ar=[]
        
        for chrom_name, chrom_data in chrom_dict.items():
            
            rep_chrom_names.append(f'{rep_name} {chrom_name}')
            chrom_name_ar.append(chrom_name)
            
    df_cor_matrix=pd.DataFrame(np.nan, index=rep_chrom_names, columns=rep_chrom_names)
    
    # Calculate track correlations and add values to the correlation matrix.
    for chrom_name in chrom_name_ar:
        
        for rep_name_1, chrom_dict_1 in Rep_chrom_dict_masked.items():
            
            track_1=chrom_dict_1[chrom_name]
            
            for rep_name_2, chrom_dict_2 in Rep_chrom_dict_masked.items():
                
                track_2=chrom_dict_2[chrom_name]
		                
                if cor_method=="pearson":
					
                    cor_coef=sci.pearsonr(track_1, track_2)[0]
                    
                elif cor_method=="spearman":
                    
                    cor_coef=sci.spearmanr(track_1, track_2)[0]
                    
                df_cor_matrix.at[f'{rep_name_1} {chrom_name}', f'{rep_name_2} {chrom_name}']=cor_coef
    
    # Save correlation matrix to file.            
    df_cor_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
                
    # Create correlation heatmap.
    fig=plt.figure(figsize=(10,10), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1, aspect="equal")
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    #Label ticks.
    labels=list(df_cor_matrix)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    #Create text annotation for heatmap pixels.
    for i in range(len(labels)):
        for j in range(len(labels)):
            text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 2), ha="center", va="center", color="black")    
    #Add colorbar.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, location='right', ticks=[-1.0, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(0.5, 1.0))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=400, figsize=(15, 10))
    plt.show()
    plt.savefig(outpath_folder+file_name+'.svg', dpi=400, figsize=(15, 10)) 
    plt.show()

    plt.close()
    
    return df_cor_matrix


#########
##Plot correlation matrix.
#########

def correlation_matrix_plot(df_cor_matrix, cor_method, title, outpath_folder, file_name):
    
    fig=plt.figure(figsize=(10,10), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    #Create correlation matrix and heatmap.
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1, aspect="equal")
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    #Label ticks.
    labels=list(Rep_chrom_dict_masked)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    #Create text annotation for heatmap pixels.
    for i in range(len(labels)):
        for j in range(len(labels)):
            text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 2), ha="center", va="center", color="black")    
    #Add colorbar.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, location='right', ticks=[-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(0.5, 1.0))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=400, figsize=(15, 10))
    plt.show()
    plt.savefig(outpath_folder+file_name+'.svg', dpi=400, figsize=(15, 10)) 
    plt.show()
    plt.close()
    
    return


#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#######

def Clustering(clust_matrix, outpath_folder, file_name):
    X = clust_matrix.values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [clust_matrix.columns.tolist()[i] for i in list((np.argsort(ind)))]
    clust_matrix = clust_matrix.reindex(columns, axis=1)
    clust_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    return clust_matrix

def read_wig_correlate_plot(dict_of_tracks, dict_of_deletions, corr_method, PWD_out):
    
    Rep_chrom_dict=read_combine(dict_of_tracks)
    Rep_chrom_dict_masked=mask_regions(Rep_chrom_dict, dict_of_deletions)
    Correlation_matrix=make_correlation_matrix_plot(Rep_chrom_dict_masked, corr_method, 'Correlation of samples', PWD_out, "TAIR10.1_tracks_correlation_matrix_masked")
    Correlation_matrix_clusterized=Clustering(Correlation_matrix, PWD_out, "TAIR10.1_tracks_correlation_matrix_clusterized_masked")
    correlation_matrix_plot(Correlation_matrix_clusterized, corr_method, 'Correlation of samples clusterized', PWD_out, "TAIR10.1_tracks_correlation_matrix_clusterized_masked")    
    
    return

read_wig_correlate_plot(Dict_of_tracks, Dict_of_deletions, Correlation_method, PWD_out)


