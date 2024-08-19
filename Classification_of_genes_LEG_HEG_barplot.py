###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis Topo-Seq analysis##

# Classification of genes: LEG, HEG.
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

# Path to genes data.
Genes_data=os.path.join(PWD, "MtRNA_Seq_BK010421.1.FPKM_combined.xlsx")

# Gene groups.
Genes_groups=["All_genes", "HEG_25", "LEG_25"]

# Output path.
Output_path=os.path.join(PWD, "BK010421.1_genes_classification")
if not os.path.isdir(Output_path):
    os.mkdir(Output_path)
    
    
    
def read_rna_seq(rna_seq_path, sheet_name):
    
    RNA_seq_df=pd.read_excel(rna_seq_path, sheet_name=sheet_name, header=0)
    
    gene_name=list(RNA_seq_df['accession'])
    FPKM_ar=list(RNA_seq_df['FPKM_average'])
    
    return gene_name, FPKM_ar


def genes_groups_barplot(RNA_Seq_data_dict, genes_groups, output_path):
    
    fig=plt.figure(figsize=(4,3), dpi=100)
    plot=fig.add_subplot(111) 
    
    All_gene_names=RNA_Seq_data_dict[genes_groups[0]][0]
    FPKM_log_ar=np.log2(np.array(RNA_Seq_data_dict[genes_groups[0]][1]))
    Coordinates_genes=range(0, len(All_gene_names), 1)
    FPKM_log_ar_sorted, All_gene_names_sorted=zip(*sorted(zip(FPKM_log_ar, All_gene_names)))
    plot.bar(Coordinates_genes, FPKM_log_ar_sorted, color="#333738")
    
    LEG_FPKM_log_ar=np.log2(np.array(RNA_Seq_data_dict[genes_groups[2]][1]))
    Coordinates_genes_LEG=range(0, len(LEG_FPKM_log_ar), 1)
    LEG_FPKM_log_ar_sorted=sorted(LEG_FPKM_log_ar)    
    plot.bar(Coordinates_genes_LEG, LEG_FPKM_log_ar_sorted, color="#757d8b")
    
    HEG_FPKM_log_ar=np.log2(np.array(RNA_Seq_data_dict[genes_groups[1]][1]))
    Coordinates_genes_HEG=range(len(All_gene_names)-len(HEG_FPKM_log_ar), len(All_gene_names), 1)
    HEG_FPKM_log_ar_sorted=sorted(HEG_FPKM_log_ar)    
    plot.bar(Coordinates_genes_HEG, HEG_FPKM_log_ar_sorted, color="#b08642")    
    
    plot.set_xlim([0-0.5,len(All_gene_names_sorted)+0.5])
    plot.set_xticks(Coordinates_genes)
    
    All_gene_names_sorted_sparsed=[]
    for i in range(len(All_gene_names_sorted)):
        if i%2==0:
            All_gene_names_sorted_sparsed.append(All_gene_names_sorted[i])
        else:
            All_gene_names_sorted_sparsed.append("")
        
    plot.set_xticklabels(All_gene_names_sorted_sparsed, rotation=90, size=5)
    plot.set_ylabel('RNA-Seq FPKM, log2', size=12)
    plot.spines["top"].set_visible(False)
    plot.spines["right"].set_visible(False)        
    plot.legend(fontsize=11, frameon=False, loc="upper right", labelspacing=0.3, handlelength=1, handletextpad=0.5)    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'Chloroplast_genes_expression_log.png'), dpi=400, figsize=(4, 3))
    plt.savefig(os.path.join(output_path, f'Chloroplast_genes_expression_log.svg'), dpi=400, figsize=(4, 3))    
    
    return


def wrapper_func(rna_seq_path, genes_groups, output_path):
    
    # Read RNA-Seq data.
    RNA_Seq_data_dict={}
    for sheet_name in genes_groups:
        gene_name_ar, FPKM_av_ar=read_rna_seq(rna_seq_path, sheet_name)
        RNA_Seq_data_dict[sheet_name]=[gene_name_ar, FPKM_av_ar]
    
    # Correlate RNA-Seq and ChIP-Seq.
    genes_groups_barplot(RNA_Seq_data_dict, genes_groups, output_path)
    
    
    return

wrapper_func(Genes_data, Genes_groups, Output_path)