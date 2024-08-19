###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase Topo-Seq analysis##

#The script takes raw GCSs data, returns only trusted GCSs, 
#computes GCSs shared between different conditions, 
#draws Venn diagrams of the sets overlappings, 
#writes GCSs sets.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import collections
from matplotlib_venn import venn2, venn3, venn3_circles, venn2_circles
import numpy as np

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Path to the working directory
pwd="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\"

#Input data
Path_to_chr1_replicas={'Rep_1_chr1': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002684.1_called.txt",
                       'Rep_2_chr1': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002684.1_called.txt",
                       }
Path_to_chr1_replicas_no_UBQ13={'Rep_1_chr1': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002684.1_called_no_UBQ13.txt",
                                'Rep_2_chr1': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002684.1_called_no_UBQ13.txt",
                                }
Path_to_chr2_replicas={'Rep_1_chr2': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002685.1_called.txt",
                       'Rep_2_chr2': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002685.1_called.txt",
                       }
Path_to_chr3_replicas={'Rep_1_chr3': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002686.1_called.txt",
                       'Rep_2_chr3': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002686.1_called.txt",
                       }
Path_to_chr4_replicas={'Rep_1_chr4': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002687.1_called.txt",
                       'Rep_2_chr4': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002687.1_called.txt",
                       }
Path_to_chr5_replicas={'Rep_1_chr5': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002688.1_called.txt",
                       'Rep_2_chr5': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002688.1_called.txt",
                       }
Path_to_chlo_replicas={'Rep_1_chlo': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_AP000423.1_called.txt",
                       'Rep_2_chlo': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_AP000423.1_called.txt",
                       }
Path_to_mito_replicas={'Rep_1_mito': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_BK010421.1_called.txt",
                       'Rep_2_mito': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_BK010421.1_called.txt",
                       }
Path_to_all_replicas={'Rep_1_chr1': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002684.1_called_no_UBQ13.txt",
                      'Rep_2_chr1': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002684.1_called_no_UBQ13.txt",
                      'Rep_1_chr2': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002685.1_called.txt",
                      'Rep_2_chr2': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002685.1_called.txt",
                      'Rep_1_chr3': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002686.1_called.txt",
                      'Rep_2_chr3': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002686.1_called.txt",
                      'Rep_1_chr4': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002687.1_called.txt",
                      'Rep_2_chr4': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002687.1_called.txt",
                      'Rep_1_chr5': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_CP002688.1_called.txt",
                      'Rep_2_chr5': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_CP002688.1_called.txt",
                      'Rep_1_chlo': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_AP000423.1_called.txt",
                      'Rep_2_chlo': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_AP000423.1_called.txt",
                      'Rep_1_mito': pwd + "GCS_calling\Gyrase_Cfx_1_GCSs_calling_0_01\Gyrase_Cfx_1_GCSs_calling_0_01Gyrase_Cfx_1_raw_dsCS_BK010421.1_called.txt",
                      'Rep_2_mito': pwd + "GCS_calling\Gyrase_Cfx_2_GCSs_calling_0_01\Gyrase_Cfx_2_GCSs_calling_0_01Gyrase_Cfx_2_raw_dsCS_BK010421.1_called.txt",
                      }


#Configuration of the output for the GCSs data in replicas.
Replicas_path_out=os.path.join(pwd, 'Trusted_GCSs_0.01')
if not os.path.exists(Replicas_path_out):
    os.makedirs(Replicas_path_out)


#######
#Parsing raw GCSs coordinates, returns dictionary - GCSs_coordinate:N3E.
#######

def read_GCSs_file(GCSs_file_path):
    GCSs_dict={}
    GCSs_in=open(GCSs_file_path, 'r')
    for line in GCSs_in:
        line=line.rstrip().split('\t')
        if line[0] not in ['dsCS_coordinate']:
            GCSs_dict[int(line[0])]=float(line[1])
    GCSs_in.close()
    return GCSs_dict

#######
#Combines replicates into one GCSs table.
#######

def combine_replicates(replicas_dict, path_out, name):
    #Merges a range of replicates
    GCSs_replicas_dict={}
    names_ar=[]
    for key, value in replicas_dict.items(): #Iterates replicas
        names_ar.append(key)
        #Read file with raw GCSs
        Raw_GCSs_dict=read_GCSs_file(value)
        for k, v in Raw_GCSs_dict.items(): #Iterates raw GCSs
            #Table filling process initiation
            if len(names_ar)==1:
                GCSs_replicas_dict[k]=[v]
            #Table filling process continuing (the table already contains at least one GCSs set)
            else:
                #If GCSs is already in the table
                if k in GCSs_replicas_dict:
                    GCSs_replicas_dict[k].append(v)
                #If this is the first occurrence of the element in a NON empty table.
                else:
                    add_el=[]
                    for j in range(len(names_ar)-1):
                        add_el.append(0)
                    add_el.append(v)
                    GCSs_replicas_dict[k]=add_el
        #If table body line contains less elements than header does, hence add zero.
        for k, v in GCSs_replicas_dict.items():
            if len(v)<len(names_ar):
                GCSs_replicas_dict[k].append(0)
    #Sorting the list of dictionary keys.
    GCSs_replicas_dict_sorted=collections.OrderedDict(sorted(GCSs_replicas_dict.items()))
    #Writes merged GCSs data
    fileout=open(os.path.join(path_out, f'{name}_GCSs_replicates.txt'), 'w')
    #Header
    fileout.write('GCSs_coordinate\t')
    for i in names_ar:
        fileout.write(str(i) + '_N3E\t')
    fileout.write('\n')
    #Body of the table
    for k, v in GCSs_replicas_dict_sorted.items():
        fileout.write(str(k) + '\t')
        for i in GCSs_replicas_dict_sorted[k]:
            fileout.write(str(i) + '\t')
        fileout.write('\n')
    fileout.close()
    
    return GCSs_replicas_dict
 

#######
#Returns only trusted GCSs - observed at least 2 times within 3 biological replicates.
#Data organization: 1. coordinate of GCSs, 2.-4. N3E values for biological replicates 1-3
#######

def trusted(ar):
    av_height=0
    ind=0
    for i in range(len(ar)):
        if ar[i]>0:
            ind=ind+1
            av_height=av_height+ar[i]
    if ind>1:
        return av_height/ind
    else:
        return "No signal"

def trusted_GCSs_calling(GCSs_dictionary):
    ar=[]
    for k, v in GCSs_dictionary.items():
        if trusted(v)!="No signal":
            ar.append([k, trusted(v)])
    return ar

def replicas_comb_trust_wrapper(replicas_dict, path_out, name):
    print('Now working with: ' + str(name))
    cur_GCSs_dict=combine_replicates(replicas_dict, path_out, name)
    cur_GCSs_trusted=trusted_GCSs_calling(cur_GCSs_dict)
    print('Number of trusted GCSs for ' + str(name) + ' : ' + str(len(cur_GCSs_trusted)))
    return cur_GCSs_trusted


#######
#Parses replicas, overlaps lists of GCSs, output data for Venn diagram construction.
#######

def replicates_parsing_to_list_and_overlapping(replicas_dict):
    #Parsing
    GCSs_dict={}
    for k, v in replicas_dict.items(): #Iterate replicas.
        GCSs_dict[k]=[]
        for c, h in read_GCSs_file(v).items(): #Iterate GCSs.
            GCSs_dict[k].append(c)
    return GCSs_dict

#######
#Venn diagram represents GCSs sets overlapping.
#description2: one, two, one_two
#description3: one, two, one_two, three, one_three, two_three, one_two_three
#######

def create_venn_2(path_to_replicas, plot_outpath, chr_set_name):

    Dict_of_replicas=replicates_parsing_to_list_and_overlapping(path_to_replicas)
    
    list_of_set_names=[]
    list_of_sets=[]
    for set_name, GCSs_rep_list in Dict_of_replicas.items():
        list_of_set_names.append(set_name)
        list_of_sets.append(set(GCSs_rep_list))
    
    plt.figure(figsize=(3.2, 1.8), dpi=100)       
    venn2(list_of_sets, set_labels=list_of_set_names)
    venn2_circles(list_of_sets)
    plt.savefig(os.path.join(plot_outpath, f'{chr_set_name}_replicas_venn.png'), dpi=320)
    plt.savefig(os.path.join(plot_outpath, f'{chr_set_name}_replicas_venn.svg'), dpi=320)
    plt.close()
    
    return

#######
#GCSs sets average N3E estimation.
#######

def average_height(ar):
    
    av_he=0
    for i in range(len(ar)):
        peak_he=np.mean(ar[i][1:])
        av_he=av_he+peak_he
        
    return av_he/len(ar)


#######
#Write down files with GCSs lists - trusted or shared.
#######

def write_trusted_GCSs(ar, output_path, chr_set_name):
    
    fileout=open(os.path.join(output_path, f'{chr_set_name}_trusted_GCSs.txt'), 'w')
    fileout.write(f'GCSs_coordinate\tRep1_N3E\tRep2_N3E\n')
    ar.sort(key=lambda tup: tup[0])

    for i in range(len(ar)):
        fileout.write(f'{ar[i][0]}\t{ar[i][1]}\n')
    fileout.close()
    
    return
    

#######
#Wrapper function.
#######

def wrapper_function(path_to_all_replicas, output_path, all_conditions_name, path_to_chr_replicas, chr_set_name):
    
    # Prepares TCSs table for all Topo-Seq conditions
    combine_replicates(path_to_all_replicas, output_path, all_conditions_name)   
    
    # Identify trusted TCSs for Topo-Seq replicates.
    GCSs_list=replicas_comb_trust_wrapper(path_to_chr_replicas, output_path, chr_set_name)
    
    # Venn diagram representing GCSs sets overlapping.
    create_venn_2(path_to_chr_replicas, output_path, chr_set_name)
    
    print(f'Average GCSs N3E for {chr_set_name}: {average_height(GCSs_list)}')

    #Write down files with GCSs lists - trusted or shared.
    write_trusted_GCSs(GCSs_list, output_path, chr_set_name)
    
    return

All_conditions_name="All_replicates_and_chromosomes_GCSs"
Chr1_set_name="CP002684.1_Chr1_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_chr1_replicas, Chr1_set_name)
Chr1_set_name_no_UBQ13="CP002684.1_no_UBQ13_Chr1_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_chr1_replicas_no_UBQ13, Chr1_set_name_no_UBQ13)
Chr2_set_name="CP002685.1_Chr2_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_chr2_replicas, Chr2_set_name)
Chr3_set_name="CP002686.1_Chr3_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_chr3_replicas, Chr3_set_name)
Chr4_set_name="CP002687.1_Chr4_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_chr4_replicas, Chr4_set_name)
Chr5_set_name="CP002688.1_Chr5_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_chr5_replicas, Chr5_set_name)
Chlo_set_name="AP000423.1_Chloroplast_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_chlo_replicas, Chlo_set_name)
Mito_set_name="BK010421.1_Mitochondria_GSCs"
wrapper_function(Path_to_all_replicas, Replicas_path_out, All_conditions_name, Path_to_mito_replicas, Mito_set_name)
 
print('Script ended its work succesfully!') 
 
