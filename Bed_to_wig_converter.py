###############################################
##Dmitry Sutormin, 2019##
##Topo-Seq/ChIP-Seq analysis##

####
#Converts bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the working directory.
PWD="/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/Mito_RNA_Seq_Garcia_2021/"
#Path to the input file
filein_path_dict={'1' :  PWD + "Cov_depth/ERR1665215.bed",  
		  '2' :  PWD + "Cov_depth/ERR1665219.bed",
		  '3' :  PWD + "Cov_depth/ERR1665220.bed",
		 }              


#Path to the output file.
fileout_path_dict={'1' :  PWD + "WIG/MtRNA_Seq_1_ERR1665215.wig",  
		   '2' :  PWD + "WIG/MtRNA_Seq_2_ERR1665219.wig",
		   '3' :  PWD + "WIG/MtRNA_Seq_3_ERR1665220.wig",
                  }

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "MtRNA_Seq_1_ERR1665215",  
	   '2' :  "MtRNA_Seq_2_ERR1665219",
           '3' :  "MtRNA_Seq_3_ERR1665220",
           }     

#ID of chromosome (for TAIR10: GCA_000001735.2_TAIR10.1)
Chromosome_name=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)


def read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual):
    for sample_name, sample_path in filein_path_dict.items():
        print(f'Now is processing: {sample_path}')
        print(f'Progress: {sample_name}/{len(filein_path_dict)}')
        
        filein=open(filein_path_dict[sample_name], 'r')
        fileout=open(fileout_path_dict[sample_name], 'w')
        
        Ar_of_Cromosome_names=[]
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in Ar_of_Cromosome_names:
                if Auto_or_manual==0:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
                elif Auto_or_manual==1:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
                Ar_of_Cromosome_names.append(line[0])
            else:
                fileout.write(line[2]+'\n')
            
        filein.close()
        fileout.close()    
    return


read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual)
                        