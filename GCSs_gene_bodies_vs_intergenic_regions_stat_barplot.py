###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase Topo-Seq analysis##

# The script analysis sets of genome intervals (genes)
# for the enrichment of GCSs (binomial test) in intergenic regions vs gene bodies.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from scipy.stats import binom

#######
#Variables to be defined.
#######

# Path to working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\Final_tracks\\"

# Input data - GCSs, TAB.
Path_to_GCSs_files={'Cfx 0.01': os.path.join(PWD, "AP000423.1_Trusted_GCSs_0.01.BroadPeak"),
                    'Cfx 0.05': os.path.join(PWD, "AP000423.1_Trusted_GCSs_0.05.BroadPeak"),
                    }

# Input data - genes, TAB.
Transcription_data_path={'HEG_50'        : os.path.join(PWD, 'Genes', 'AP000423.1_no_exon_HETU_50.bed'),
                         'All_genes_133'  : os.path.join(PWD, 'Genes', 'AP000423.1_no_exon.bed'),
                         'LEG_50'        : os.path.join(PWD, 'Genes', 'AP000423.1_no_exon_LETU_50.bed'),
                         }

# Genome length, bp.
Genome_len=154478     # A. thaliana chloroplast genome: 154478 bp, A. thaliana mitochondrial genome: 367808 bp

# Output path.
Output_path=os.path.join(PWD, "Gene_body_IGR_analysis", "AP000423.1")
if not os.path.isdir(Output_path):
    os.mkdir(Output_path) 
    
    
#######
#Trusted GCSs data parsing.
#######

def trusted_GCSs_parsing(input_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        ar=[]
        filein=open(v, 'r')
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in ['GCSs_coordinate']:
                ar.append(int(line[1]))
            else:
                continue
        GCSs_sets_dict[k]=ar
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(ar)))
        
    return GCSs_sets_dict  


#######
#Parsing genes data.
#######

def TUs_parser(TUs_sets_path):
    
    TUs_sets={}
    TUs_mean_len_dict={}
    for k, v in TUs_sets_path.items():
        filein=open(v, 'r')
        plus=0
        minus=0
        ar=[]
        TU_len_ar=[]
        for line in filein:
            tu={}
            line=line.rstrip().split('\t')
            if line[0] not in ["GeneID", "OperonID"]:
                tu['TUID']=str(line[3]) #TUID
                tu['TU name']=str(line[3]) #TU name/composition
                tu['Start']=int(line[1]) #Start
                tu['End']=int(line[2]) #End
                tu['Strand']=str(line[5]) #Strand +/-
                tu['Transcription level']=float(line[4].replace(',','.')) #Transcription level
                ar.append(tu)
                TU_len_ar.append(int(line[2])-int(line[1]))
                if line[5]=='+':
                    plus+=1
                elif line[5]=="-":
                    minus+=1
        filein.close()
        
        TUs_sets[k]=ar
        TUs_mean_len_dict[k]=np.mean(TU_len_ar)
        
        print(f'Number of TUs in forward for {k} set: {plus}')
        print(f'Number of TUs in reverse for {k} set: {minus}')
        print(f'Mean length of TUs {k} set: {np.mean(TU_len_ar)}')
        
    return TUs_sets, TUs_mean_len_dict


def TU_IGR_association(GSCs_data_dict, gene_set_data, output_path, gene_set_name, genome_len):
    
    GCSs_genes_and_IGR_assoc_norm_kb={}
    
    Genome_track=[0]*genome_len
    
    for gene_info in gene_set_data:
        
        Gene_start=gene_info['Start']
        Gene_end=gene_info['End']
        
        Genome_track[Gene_start:Gene_end]=[1]*(Gene_end-Gene_start)
        
    Sum_gene_len=sum(Genome_track)
    Sum_IGR_len=genome_len-Sum_gene_len
    GCSs_genes_and_IGR_assoc_norm_kb["Genes_len"]=Sum_gene_len
    GCSs_genes_and_IGR_assoc_norm_kb["IGR_len"]=Sum_IGR_len
    
    print(f'Total length of gene bodies for {gene_set_name}: {Sum_gene_len} bp')
    print(f'Total length of IGRs for {gene_set_name}: {Sum_IGR_len} bp')
    
    for GCSs_set, GSCs_ar in GSCs_data_dict.items():
        
        Gene_GCSs_num=0
        IGR_GCSs_num=0
        
        for GCSs_coord in GSCs_ar:
            
            if sum(Genome_track[GCSs_coord:GCSs_coord+4])>0:
                Gene_GCSs_num+=1
            else:
                IGR_GCSs_num+=1
                
        print(f'Number of GCSs {GCSs_set} observed in gene bodies for {gene_set_name}: {Gene_GCSs_num}')
        print(f'Number of GCSs {GCSs_set} observed in IGRs for {gene_set_name}: {IGR_GCSs_num}')
        
        Gene_GCSs_num_per_kb=Gene_GCSs_num/(Sum_gene_len/1000)
        IGR_GCSs_num_per_kb=IGR_GCSs_num/(Sum_IGR_len/1000)
        
        print(f'Number of GCSs {GCSs_set} per kb observed in gene bodies for {gene_set_name}: {Gene_GCSs_num_per_kb}/kb')
        print(f'Number of GCSs {GCSs_set} per kb observed in IGRs for {gene_set_name}: {IGR_GCSs_num_per_kb}/kb')  
        
        print(f'Statistics for GCSs {GCSs_set} observed in gene bodies for {gene_set_name}: {binom.cdf(Gene_GCSs_num, len(GSCs_ar), Sum_gene_len/genome_len)}')
        print(f'Statistics for GCSs {GCSs_set} observed in IGRs for {gene_set_name}: {binom.cdf(IGR_GCSs_num, len(GSCs_ar), Sum_IGR_len/genome_len)}')
        
        GCSs_genes_and_IGR_assoc_norm_kb[GCSs_set]={"Gene_GCS_num" : Gene_GCSs_num, 
                                                    "IGR_GCS_num" : IGR_GCSs_num,
                                                    "Gene_GCS_num_per_kb" : Gene_GCSs_num_per_kb, 
                                                    "IGR_GCS_num_per_kb" : IGR_GCSs_num_per_kb,}
        
    return GCSs_genes_and_IGR_assoc_norm_kb


#######
#Wrapper function.
#######

def wrapper_func(path_to_GCSs_files, transcription_data_path, genome_len, output_path):
    
    # Read GCSs data.
    GSCs_data_dict=trusted_GCSs_parsing(path_to_GCSs_files)
    
    # Read TUs data.
    TUs_data_dict, TUs_mean_len_dict=TUs_parser(transcription_data_path)
    
    # Genes-association analysis.
    GCSs_num_norm_kb_dict={}
    for gene_set_name, gene_set_data in TUs_data_dict.items():
        GCSs_genes_and_IGR_assoc_norm_kb=TU_IGR_association(GSCs_data_dict, gene_set_data, output_path, gene_set_name, genome_len)
        GCSs_num_norm_kb_dict[gene_set_name]=GCSs_genes_and_IGR_assoc_norm_kb
        
        #TU_interval_stat_analysis(GSCs_data_dict, GCSs_all_genes_assoc_info, gene_set_data, window_width, 'genes', output_path, gene_set_name, genome_len)
        #GCSs_set_exp_interval_dict_ag=GCSs_number_norm(GCSs_all_genes_assoc_info, GSCs_data_dict, genome_len)
        #write_GCSs_norm(GCSs_set_exp_interval_dict_ag, output_path, gene_set_name)    
    
    # Plot normalized numbers of GCSs associated with US, GB, DS regions.
    #plot_GCSs_numbers(GCSs_num_norm_kb_dict, output_path)
    
    return
    
wrapper_func(Path_to_GCSs_files, Transcription_data_path, Genome_len, Output_path)

