###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase Topo-Seq analysis##

# The script analysis sets of genome intervals (genes)
# for the enrichment of GCSs (binomial test).
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
Path_to_GCSs_files={'Cfx 0.01': os.path.join(PWD, "BK010421.1_Trusted_GCSs_0.01.BroadPeak"),
                    'Cfx 0.05': os.path.join(PWD, "BK010421.1_Trusted_GCSs_0.05.BroadPeak"),
                    }

# Input data - genes, TAB.
Transcription_data_path={'HEG_25'        : os.path.join(PWD, 'Genes', 'BK010421.1_no_exon_HETU_25.bed'),
                         'All_genes_61'  : os.path.join(PWD, 'Genes', 'BK010421.1_no_exon.bed'),
                         'LEG_25'        : os.path.join(PWD, 'Genes', 'BK010421.1_no_exon_LETU_25.bed'),
                         }

# Genome length, bp.
Genome_len=154478

# Output path.
Output_path=os.path.join(PWD, "US_GB_DS_analysis", "BK010421.1")
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

#######
#GCSs association with TUs.
#######

def TU_association(GCSs_sets_dict, TU_set, set_type, window_width, path_out, set_name):
    
    fileout=open(os.path.join(path_out, f'{set_name}_numbers_of_associated_GCSs.txt'), 'w')
    if set_type=='genes':
        fileout.write("Gene_ID\tGenes_ID\tStart\tEnd\tStrand\tExpression\t") 
        for k in range(len(GCSs_sets_dict)):
            fileout.write("Condition\tGCSs in US\tGCSs in GB\tGCSs in DS\t")
        fileout.write("\n")
    else:
        print('Unknown type of the set! Possible types: 16S_operons, operons, genes.')
        return
    
    operons_stat_norm_kb={}
    operons_stat_ds={} #{condition : [number of GCSs, mean N3E, mean score, len of regions]}
    ds_regions_len=0
    
    for j in TU_set: #j - particular operon info (dictionary)  
        fileout.write(j['TUID'] + '\t' + j['TU name'] + '\t' + str(j['Start']) + '\t' + str(j['End']) + '\t' + j['Strand'] + '\t' + str(j['Transcription level']) + '\t')      
        ds_regions_len+=window_width
        for a, s in GCSs_sets_dict.items(): #a - Topo-Seq condition, s - corresponding set of GCSs.
            if a not in operons_stat_ds:
                if set_type=='genes':
                    operons_stat_ds[a]=[0, [], [], 0, [], [], 0, [], []] #[number of GCSs, N3E values, score values] for US, GB, DS correspondingly.
            fileout.write(a + '\t')
            stats=[0, 0, 0, 0, 0] #USUS USGB GB GBDS DSDS

            for k in s: 
                if j['Start']+window_width>k>j['Start'] and j['Strand']=='-': #GBDS
                    stats[3]+=1
                  
                if j['Start']+window_width>k>j['Start'] and j['Strand']=='+': #USGB
                    stats[1]+=1

                if j['End']>k>j['End']-window_width and j['Strand']=='-': #USGB
                    stats[1]+=1
  
                if j['End']>k>j['End']-window_width and j['Strand']=='+': #GBDS
                    stats[3]+=1
                
                if j['Start']>k>j['Start']-window_width and j['Strand']=='-': #DSDS
                    stats[4]+=1

                if j['Start']>k>j['Start']-window_width and j['Strand']=='+': #USUS
                    stats[0]+=1   

                if j['End']+window_width>k>j['End'] and j['Strand']=='-': #USUS
                    stats[0]+=1

                if j['End']+window_width>k>j['End'] and j['Strand']=='+': #DSDS
                    stats[4]+=1 

                if j['End']>k>j['Start']: #GB
                    stats[2]+=1

            if set_type=='genes':
                fileout.write(str(stats[0]) + '\t' + str(stats[2]) + '\t' + str(stats[4]) + '\t') #US, GB, DS
                #US
                operons_stat_ds[a][0]+=stats[0]
                #GB
                operons_stat_ds[a][3]+=stats[2]
                #DS
                operons_stat_ds[a][6]+=stats[4]
                
            US_norm=stats[0]/(window_width/1000)
            GB_norm=stats[2]/((j['End']-j['Start'])/1000)
            DS_norm=stats[4]/(window_width/1000)
            
            if a not in operons_stat_norm_kb:
                operons_stat_norm_kb[a]=[[US_norm], [GB_norm], [DS_norm]]
            else:
                operons_stat_norm_kb[a][0].append(US_norm)
                operons_stat_norm_kb[a][1].append(GB_norm)
                operons_stat_norm_kb[a][2].append(DS_norm)
            
                 
        fileout.write('\n')
    for a, v in operons_stat_ds.items():
        v.append(ds_regions_len)
    fileout.close()
    
    return operons_stat_ds, operons_stat_norm_kb #{condition : [##number of GCSs, N3E values, score values##*1-4, len of regions]}


#######
#Statistical analysis of intevals:
#1) Number of GCSs under/overrepresentation - binomail test.
#2) Interval GCSs mean N3E value vs all GCSs mean N3E - t-test.
#3) Interval GCSs mean score value vs all GCSs mean score - t-test.
#4) Interval mean score vs genome mean score - t-test.
#######

def TU_interval_stat_analysis(GCSs_sets_dict, intervals_GCSs_dict, intervals, window_width, set_type, path_out, set_name, genome_len):
    
    #Correcting genome length.
    deletions=[] #Deletions.
    del_len=0
    for i in deletions:
        del_len+=i[1]-i[0]
    genome_len_dc=genome_len-del_len
            
    fileout=open(os.path.join(path_out, f'{set_name}_compartments_statistics_GCSs_number.txt'), 'w')
    #Performes t-test for comparison of intervals GCSs N3E and scores with all GCSs N3E and scores.
    #Performes binomail test for enrichment estimation of GCSs fall into intervals.
    fileout.write('Test\tAntibiotic\t')
    for i in range(len(intervals_GCSs_dict.values())-1):
        fileout.write('Examined value\tValue (intervals)\tValue (overall)\tp-value\tAdditional\t')
    fileout.write('\n') 
    
    if set_type=='genes':
        for a, ns in GCSs_sets_dict.items():

            #Number of GCSs
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][0], len(ns), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in US
            fileout.write('binomial test\t' + a + '\tUS Number of GCSs\t' + str(intervals_GCSs_dict[a][0]) + '\t' + 
                          str(len(ns)) + '\t' + str(GCSs_number_stat) + '\n')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][3], len(ns), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in GB
            fileout.write('binomial test\t' + a + '\tGB Number of GCSs\t' + str(intervals_GCSs_dict[a][3]) + '\t' + 
                          str(len(ns)) + '\t' + str(GCSs_number_stat) + '\n')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][6], len(ns), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in DS
            fileout.write('binomial test\t' + a + '\tDS Number of GCSs\t' + str(intervals_GCSs_dict[a][6]) + '\t' + 
                          str(len(ns)) + '\t' + str(GCSs_number_stat) + '\n')            
              
    fileout.close() 
    
    return


#######
#Statistical analysis for the number of GCSs associated with genes, additional normalization step is added. 
#######

def GCSs_number_norm(intervals_GCSs_dict, GCSs_sets_dict, genome_len):
    #Correcting genome length.
    deletions=[] #Deletions.
    del_len=0
    for i in deletions:
        del_len+=i[1]-i[0]
    genome_len_dc=genome_len-del_len  
    #Number of GCSs normalization and statistics.
    GCSs_set_exp_interval_dict={}
    for a, ns in intervals_GCSs_dict.items():
        Num_GCSs_expected=len(GCSs_sets_dict[a])*ns[-1]/genome_len_dc
        GCSs_set_exp_interval_dict[a]=[Num_GCSs_expected] 
        for j in range(int((len(ns)-1)/3)):
            GCSs_set_exp_interval_dict[a].append(ns[j*3]) #Number of GCSs fall into particular compartment (US, GB, DS) of 16S operons.
            GCSs_set_exp_interval_dict[a].append(binom.cdf(ns[j*3], len(GCSs_sets_dict[a]), ns[-1]/genome_len_dc)) #p-value of binomial test for the number of GCSs fall into particular compartment (US, GB, DS) of the 16S operons.
            GCSs_set_exp_interval_dict[a].append(float(ns[j*3])/Num_GCSs_expected) #Normalized number of GCSs fall into particular compartment (US, GB, DS) of the 16S operons.
    
    return GCSs_set_exp_interval_dict #[GCSs expected] + [GCSs obs, p-value, GCSs norm]*[US, GB, DS]

  
#######
#Writes GCSs numbers information to file: Number of GCSs expected, Number of GCSs observed, p-value - binomial test, Number of GCSs normalized.
#######

def write_GCSs_norm(GCSs_set_exp_interval_dict, path_out, set_name):
    fileout=open(os.path.join(path_out, f'{set_name}_normalized_GCSs_numbers_and_statistics.txt'), 'w')
    fileout.write('Condition\tCompartment\tNumber of GCSs expected\tNumber of GCSs observed\tp-value\tNumber of GCSs normalized\n')
    Compartment_names=['US', 'GB', 'DS']

    for a, s in GCSs_set_exp_interval_dict.items():
        for i in range(len(Compartment_names)):
            fileout.write(a + '\t' + Compartment_names[i] + '\t' + str(round(s[0],3)) + '\t' + str(s[(i*3)+1]) + '\t' + str(s[(i*3)+2]) + '\t' + str(round(s[(i*3)+3],3)) +'\n')
    fileout.close()
    return


def simple_beeswarm(y, nbins=None): # Taken from https://stackoverflow.com/questions/36153410/how-to-create-a-swarm-plot-with-matplotlib
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """
    y = np.asarray(y)
    if nbins is None:
        nbins = len(y) // 6

    # Get upper bounds of bins
    x = np.zeros(len(y))
    ylo = np.min(y)
    yhi = np.max(y)
    dy = (yhi - ylo) / nbins
    ybins = np.linspace(ylo + dy, yhi - dy, nbins - 1)

    # Divide indices into bins
    i = np.arange(len(y))
    ibs = [0] * nbins
    ybs = [0] * nbins
    nmax = 0
    for j, ybin in enumerate(ybins):
        f = y <= ybin
        ibs[j], ybs[j] = i[f], y[f]
        nmax = max(nmax, len(ibs[j]))
        f = ~f
        i, y = i[f], y[f]
    ibs[-1], ybs[-1] = i, y
    nmax = max(nmax, len(ibs[-1]))

    # Assign x indices
    dx = 1 / (nmax // 2)
    for i, y in zip(ibs, ybs):
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(y)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x


def plot_GCSs_numbers(GCSs_num_norm_kb_dict, output_path):
    
    Topo_Seq_conditions_list=list(GCSs_num_norm_kb_dict[list(GCSs_num_norm_kb_dict.keys())[0]].keys())
    Genes_sets_list=list(GCSs_num_norm_kb_dict.keys())
    
    segments_ar=['US', 'GB', 'DS']
    color_ar=['#77d96a', '#3ec1db', '#de7976']
    coordinates=[1,2,3]
    
    coords_ar=[]
    for i in range(len(Genes_sets_list)):
        for coord in coordinates:
            coords_ar.append((i*4)+coord)
    
    for Topo_condition in Topo_Seq_conditions_list:
        
        fig=plt.figure(figsize=(4,3), dpi=100)
        plot=fig.add_subplot(111) 
        
        i=0
        labels_ar=[]
        
        for gene_set, GCSs_data_dict in GCSs_num_norm_kb_dict.items():
            
            gene_set_data=GCSs_data_dict[Topo_condition]
            
            print(Topo_condition, gene_set)
            print(np.mean(gene_set_data[0]))
            print(np.mean(gene_set_data[1]))
            print(np.mean(gene_set_data[2]))
            
            if i!=2:
                plot.bar(coords_ar[(i*3)+0], np.mean(gene_set_data[0]), color=color_ar[0], edgecolor='k')
                plot.errorbar(coords_ar[(i*3)+0], np.mean(gene_set_data[0]), yerr=np.std(gene_set_data[0]), capsize=4, elinewidth=1, fmt="none", color="k")                
                plot.bar(coords_ar[(i*3)+1], np.mean(gene_set_data[1]), color=color_ar[1], edgecolor='k')
                plot.errorbar(coords_ar[(i*3)+1], np.mean(gene_set_data[1]), yerr=np.std(gene_set_data[1]), capsize=4, elinewidth=1, fmt="none", color="k")
                plot.bar(coords_ar[(i*3)+2], np.mean(gene_set_data[2]), color=color_ar[2], edgecolor='k')
                plot.errorbar(coords_ar[(i*3)+2], np.mean(gene_set_data[2]), yerr=np.std(gene_set_data[2]), capsize=4, elinewidth=1, fmt="none", color="k")
            else:
                plot.bar(coords_ar[(i*3)+0], np.mean(gene_set_data[0]), color=color_ar[0], edgecolor='k', label='Upstream')
                plot.errorbar(coords_ar[(i*3)+0], np.mean(gene_set_data[0]), yerr=np.std(gene_set_data[0]), capsize=4, elinewidth=1, fmt="none", color="k")
                plot.bar(coords_ar[(i*3)+1], np.mean(gene_set_data[1]), color=color_ar[1], edgecolor='k', label='Gene body')
                plot.errorbar(coords_ar[(i*3)+1], np.mean(gene_set_data[1]), yerr=np.std(gene_set_data[1]), capsize=4, elinewidth=1, fmt="none", color="k")
                plot.bar(coords_ar[(i*3)+2], np.mean(gene_set_data[2]), color=color_ar[2], edgecolor='k', label='Downstream')
                plot.errorbar(coords_ar[(i*3)+2], np.mean(gene_set_data[2]), yerr=np.std(gene_set_data[2]), capsize=4, elinewidth=1, fmt="none", color="k")
            
            i+=1
            labels_ar+=['', gene_set,]

        plot.set_xticklabels(labels_ar)
        plot.set_ylabel('GCSs/kb', size=12)
        plot.spines["top"].set_visible(False)
        plot.spines["right"].set_visible(False)        
        plot.legend(fontsize=11, frameon=False, loc="upper right", labelspacing=0.3, handlelength=1, handletextpad=0.5)    
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f'GCSs_norm_{Topo_condition.replace(" ", "_")}.png'), dpi=400, figsize=(4, 3))
        plt.savefig(os.path.join(output_path, f'GCSs_norm_{Topo_condition.replace(" ", "_")}.svg'), dpi=400, figsize=(4, 3))  
        
    for Topo_condition in Topo_Seq_conditions_list:
        
        for gene_set_1, GCSs_data_dict_1 in GCSs_num_norm_kb_dict.items():
            
            gene_set_data_1=GCSs_data_dict_1[Topo_condition]
            
            for i in range(len(segments_ar)):
                
                set_1=gene_set_data_1[i]
                
                for gene_set_2, GCSs_data_dict_2 in GCSs_num_norm_kb_dict.items():
                    
                    gene_set_data_2=GCSs_data_dict_2[Topo_condition]
                    
                    for j in range(len(segments_ar)):
                        
                        set_2=gene_set_data_2[j]                
              
                        ttest_res=stats.ttest_ind(set_1, set_2)
                            
                        print(f'{Topo_condition} {gene_set_1} {segments_ar[i]} vs {gene_set_2} {segments_ar[j]}: {ttest_res}')
                        
    return
    
    
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
        window_width=TUs_mean_len_dict[gene_set_name]
        GCSs_all_genes_assoc_info, GCSs_all_genes_assoc_norm_kb=TU_association(GSCs_data_dict, gene_set_data, 'genes', window_width, output_path, gene_set_name)
        GCSs_num_norm_kb_dict[gene_set_name]=GCSs_all_genes_assoc_norm_kb
        TU_interval_stat_analysis(GSCs_data_dict, GCSs_all_genes_assoc_info, gene_set_data, window_width, 'genes', output_path, gene_set_name, genome_len)
        GCSs_set_exp_interval_dict_ag=GCSs_number_norm(GCSs_all_genes_assoc_info, GSCs_data_dict, genome_len)
        write_GCSs_norm(GCSs_set_exp_interval_dict_ag, output_path, gene_set_name)    
    
    # Plot normalized numbers of GCSs associated with US, GB, DS regions.
    plot_GCSs_numbers(GCSs_num_norm_kb_dict, output_path)
    
    return
    
wrapper_func(Path_to_GCSs_files, Transcription_data_path, Genome_len, Output_path)
