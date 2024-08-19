###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis Topo-Seq analysis##

# Correlate RNA-Seq and RNAP ChIP-Seq data.
###############################################

#######
#Packages to be imported.
#######

import os
from Bio import SeqIO
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import matplotlib.pyplot as plt

# Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\Final_tracks\\"

# Path to wig track.
Wig_track_path=os.path.join(PWD, "BK010421.1_Topo_Seq_FE_av.wig")

# Wig track name.
Wig_track_name="Gyrase FE"

# Path to RNA-Seq data.
RNA_seq_path=os.path.join(PWD, "MtRNA_Seq_BK010421.1.FPKM_combined.xlsx")

# Excel sheet name.
Sheet_name="All_genes"

# Column name.
Column_name="FPKM_average"

# Output path.
Output_path=os.path.join(PWD, "Test_BK010421.1_RNA_Seq_Topo_Seq_correlation")
if not os.path.isdir(Output_path):
    os.mkdir(Output_path)


#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    
    wigin=open(wigfile, 'r')
    NE_values_dict={}
    NE_values=[]
    
    for line in wigin:
        line=line.rstrip().split(' ')
        
        if line[0]!='track':
        
            if line[0]=='fixedStep':
                
                if len(NE_values)>0:
                    
                    NE_values_dict[chr_id]=NE_values
                
                chr_id=line[1].split('=')[1]
                NE_values=[]               
                
            else:
                
                NE_values.append(float(line[0]))
            
    NE_values_dict[chr_id]=NE_values
    
    wigin.close()
    
    return NE_values_dict


def read_rna_seq(rna_seq_path, sheet_name, column_name):
    
    RNA_seq_df=pd.read_excel(rna_seq_path, sheet_name=sheet_name, header=0)
    
    gene_name=list(RNA_seq_df['accession'])
    start_coord=list(RNA_seq_df['st'])
    end_coord=list(RNA_seq_df['end'])
    FPKM_av=list(RNA_seq_df['FPKM_average'])
    product_type=list(RNA_seq_df['Product'])
    if 'Polymerse' in RNA_seq_df.columns:
        RNAP_type=list(RNA_seq_df['Polymerse'])
    else:
        RNAP_type=['-']*len(product_type)
    
    return gene_name, start_coord, end_coord, FPKM_av, product_type, RNAP_type


def FE_for_genes(Wig_data_dict, start_coord_ar, end_coord_ar):
    
    for chrom_id, wig_track in Wig_data_dict.items():
        
        FE_mean_ar=[]
        
        for i in range(len(start_coord_ar)):
            start=start_coord_ar[i]
            end=end_coord_ar[i]
            FE_mean=np.mean(wig_track[start:end])
            
            FE_mean_ar.append(FE_mean)
    
    return FE_mean_ar


#######
#Correlation of EcTopoI and EcTopoI fold enrichments for TUs.
#######

def factors_association_EcTopoI_RpoC(FPKM_av_ar, FE_mean_ar, wig_track_name, gene_name_ar, product_type_ar, RNAP_type_ar, output_path):
    
    # No data transformation.
    factor_1=FPKM_av_ar
    factor_2=FE_mean_ar
    
    ##Linnear fitting of log data.
    fit=np.polyfit(factor_1, factor_2, 1)
    fit_fn=np.poly1d(fit) 
    xdata=np.linspace(min(factor_1), max(factor_1), 50)
    
    ##Correlation.
    pearson_cor=scipy.stats.pearsonr(factor_1, factor_2)
    print(f'Pearson correlation (RNA-Seq FPKM, {wig_track_name}) for genes: {pearson_cor}') 
    spearman_cor=scipy.stats.spearmanr(factor_1, factor_2)
    print(f'Spearman correlation (RNA-Seq FPKM, {wig_track_name}) for genes: {spearman_cor}') 
    
    ##Plot data.
    fig, plot=plt.subplots(1,1,figsize=(3,2.7), dpi=100)
    for i in range(len(product_type_ar)):
        if product_type_ar[i]=='tRNA':
            plot.scatter(factor_1[i], factor_2[i], s=2, color='green')
        elif product_type_ar[i]=='rRNA':
            plot.scatter(factor_1[i], factor_2[i], s=2, color='red')
        else:
            if RNAP_type_ar[i]=='NEP':
                plot.scatter(factor_1[i], factor_2[i], s=6, color='blue', marker='^')
            else:
                plot.scatter(factor_1[i], factor_2[i], s=2, color='blue')
       
    plot.plot(xdata, fit_fn(xdata), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot.annotate(f'Spearman cor coef: {np.round(spearman_cor[0],3)}\np-value: {"{:.1e}".format(spearman_cor[1])}', xy=(-0.5, 4), xycoords='data', size=7)
    
    for i in range(len(factor_1)):
        plot.annotate(gene_name_ar[i], xy=(factor_1[i], factor_2[i]), xytext=(factor_1[i]+0.05, factor_2[i]+0.05), xycoords='data', size=3)
    
    plot.set_xlabel('RNA-Seq FPKM', size=12)
    plot.set_ylabel(wig_track_name, size=12)
    plot.spines["top"].set_visible(False)
    plot.spines["right"].set_visible(False)
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper left')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(os.path.join(output_path, f'RNA_Seq_FPKM_vs_{wig_track_name.replace(" ", "_")}_for_genes.png'), size=(3,2.7), dpi=300) 
    plt.savefig(os.path.join(output_path, f'RNA_Seq_FPKM_vs_{wig_track_name.replace(" ", "_")}_for_genes.svg'), size=(3,2.7), dpi=300)
    
    
    # Semi-log data transformation.
    factor_1=FPKM_av_ar
    factor_2=FE_mean_ar
    
    ##Factor transformation.
    factor_1=np.log2(factor_1)    
    
    ##Linnear fitting of log data.
    fit=np.polyfit(factor_1, factor_2, 1)
    fit_fn=np.poly1d(fit) 
    xdata=np.linspace(min(factor_1), max(factor_1), 50)
    
    ##Correlation.
    pearson_cor=scipy.stats.pearsonr(factor_1, factor_2)
    print(f'Pearson correlation (log RNA-Seq FPKM, {wig_track_name}) for genes: {pearson_cor}') 
    spearman_cor=scipy.stats.spearmanr(factor_1, factor_2)
    print(f'Spearman correlation (log RNA-Seq FPKM, {wig_track_name}) for genes: {spearman_cor}') 
    
    ##Plot data.
    fig, plot=plt.subplots(1,1,figsize=(3,2.7), dpi=100)
    for i in range(len(product_type_ar)):
        if product_type_ar[i]=='tRNA':
            plot.scatter(factor_1[i], factor_2[i], s=2, color='green')
        elif product_type_ar[i]=='rRNA':
            plot.scatter(factor_1[i], factor_2[i], s=2, color='red')
        else:
            if RNAP_type_ar[i]=='NEP':
                plot.scatter(factor_1[i], factor_2[i], s=6, color='blue', marker='^')
            else:
                plot.scatter(factor_1[i], factor_2[i], s=2, color='blue')        
               
    plot.plot(xdata, fit_fn(xdata), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot.annotate(f'Spearman cor coef: {np.round(spearman_cor[0],3)}\np-value: {"{:.1e}".format(spearman_cor[1])}', xy=(-0.5, 4), xycoords='data', size=7)
    
    for i in range(len(factor_1)):
        plot.annotate(gene_name_ar[i], xy=(factor_1[i], factor_2[i]), xytext=(factor_1[i]+0.05, factor_2[i]+0.05), xycoords='data', size=3)
    
    plot.set_xlabel('RNA-Seq FPKM, log2', size=12)
    plot.set_ylabel(wig_track_name, size=12)
    plot.spines["top"].set_visible(False)
    plot.spines["right"].set_visible(False)
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper left')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(os.path.join(output_path, f'RNA_Seq_FPKM_vs_{wig_track_name.replace(" ", "_")}_for_genes_semi_log.png'), size=(3,2.7), dpi=300) 
    plt.savefig(os.path.join(output_path, f'RNA_Seq_FPKM_vs_{wig_track_name.replace(" ", "_")}_for_genes_semi_log.svg'), size=(3,2.7), dpi=300)    
    
    
    # log-log data transformation.
    factor_1=FPKM_av_ar
    factor_2=FE_mean_ar    
    
    ##Factor transformation.
    factor_1=np.log2(factor_1)
    factor_2=np.log2(factor_2)
    
    ##Linnear fitting of log data.
    fit=np.polyfit(factor_1, factor_2, 1)
    fit_fn=np.poly1d(fit) 
    xdata=np.linspace(min(factor_1), max(factor_1), 50)
    
    ##Correlation.
    pearson_cor=scipy.stats.pearsonr(factor_1, factor_2)
    print(f'Pearson correlation (log RNA-Seq FPKM, log {wig_track_name}) for genes: {pearson_cor}') 
    spearman_cor=scipy.stats.spearmanr(factor_1, factor_2)
    print(f'Spearman correlation (log RNA-Seq FPKM, log {wig_track_name}) for genes: {spearman_cor}') 
    
    ##Plot data.
    fig, plot=plt.subplots(1,1,figsize=(3,2.7), dpi=100)
    for i in range(len(product_type_ar)):
        if product_type_ar[i]=='tRNA':
            plot.scatter(factor_1[i], factor_2[i], s=2, color='green')
        elif product_type_ar[i]=='rRNA':
            plot.scatter(factor_1[i], factor_2[i], s=2, color='red')
        else:
            if RNAP_type_ar[i]=='NEP':
                plot.scatter(factor_1[i], factor_2[i], s=6, color='blue', marker='^')
            else:
                plot.scatter(factor_1[i], factor_2[i], s=2, color='blue')        
               
    plot.plot(xdata, fit_fn(xdata), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot.annotate(f'Spearman cor coef: {np.round(spearman_cor[0],3)}\np-value: {"{:.1e}".format(spearman_cor[1])}', xy=(-0.5, 4), xycoords='data', size=7)
    
    for i in range(len(factor_1)):
        plot.annotate(gene_name_ar[i], xy=(factor_1[i], factor_2[i]), xytext=(factor_1[i]+0.05, factor_2[i]+0.05), xycoords='data', size=3)
    
    plot.set_xlabel('RNA-Seq FPKM, log2', size=12)
    plot.set_ylabel(f'{wig_track_name}, log2', size=12)
    plot.spines["top"].set_visible(False)
    plot.spines["right"].set_visible(False)
    plt.legend(fontsize=8, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, loc='upper left')
    plt.tight_layout()
    plt.show()
    
    plt.savefig(os.path.join(output_path, f'RNA_Seq_FPKM_vs_{wig_track_name.replace(" ", "_")}_for_genes_log.png'), size=(3,2.7), dpi=300) 
    plt.savefig(os.path.join(output_path, f'RNA_Seq_FPKM_vs_{wig_track_name.replace(" ", "_")}_for_genes_log.svg'), size=(3,2.7), dpi=300)
    
    return



def wrapper_func(wig_track_path, wig_track_name, rna_seq_path, sheet_name, column_name, output_path):
    
    # Read wig data.
    Wig_data_dict=wig_parsing(wig_track_path)
    
    # Read RNA-Seq data.
    gene_name_ar, start_coord_ar, end_coord_ar, FPKM_av_ar, product_type_ar, RNAP_type_ar=read_rna_seq(rna_seq_path, sheet_name, column_name)
    
    # Calculate ChIP-Seq for genes.
    FE_mean_ar=FE_for_genes(Wig_data_dict, start_coord_ar, end_coord_ar)
    
    # Correlate RNA-Seq and ChIP-Seq.
    factors_association_EcTopoI_RpoC(FPKM_av_ar, FE_mean_ar, wig_track_name, gene_name_ar, product_type_ar, RNAP_type_ar, output_path)
    
    
    return

wrapper_func(Wig_track_path, Wig_track_name, RNA_seq_path, Sheet_name, Column_name, Output_path)