###############################################
##Dmitry Sutormin, 2024##
##Arabidopsis gyrase Topo-Seq analysis##

####
#Script creates and plots a web-logo for a multiple alignment of DNA sequences.
####

###############################################

#######
#Packages to be imported.
#######

import weblogo
from weblogo import LogoOptions, LogoFormat, eps_formatter, read_seq_data, LogoData, png_formatter, pdf_formatter


#######
#Creates motif logo.
#######

def Create_logo(alig_inpath, out_path):
    MFA_data=open(alig_inpath)
    MFA_seqs=read_seq_data(MFA_data)
    logodata=LogoData.from_seqs(MFA_seqs)
    logooptions=LogoOptions(yaxis_scale=1.8, pad_right=True, stacks_per_line=20)
    logooptions.show_errorbars=False
    logoformat=LogoFormat(logodata, logooptions)
    #pdf=weblogo.logo_formatter.pdf_formatter(logodata, logoformat)
    pdf=weblogo.logo_formatter.png_formatter(logodata, logoformat)
    logout=open(out_path, 'wb')
    logout.write(pdf)
    logout.close()
    return


Alig_inpath1="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\Trusted_GCSs\Motif_0.01\AP000423.1_sequences_under_GCSs_20bp_LOGO.fasta" 
Out_path1="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\Trusted_GCSs\Motif_0.01\AP000423.1_sequences_under_GCSs_20bp_LOGO.png"
Create_logo(Alig_inpath1, Out_path1)

Alig_inpath2="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\Trusted_GCSs\Motif_0.01\BK010421.1_sequences_under_GCSs_20bp_LOGO.fasta" 
Out_path2="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\NUMT_masked\\Trusted_GCSs\Motif_0.01\BK010421.1_sequences_under_GCSs_20bp_LOGO.png"
Create_logo(Alig_inpath2, Out_path2)