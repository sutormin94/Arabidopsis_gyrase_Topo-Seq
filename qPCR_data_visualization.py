###############################################
##Dmitry Sutormin, 2024##
##qPCR data visualization##

#Takes table with Ct data and primers efficiency and makes barpolts.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom



#################
### Primers calibration data analysis.
#################


#Path to qPCR data table.
plant_1_pc_table="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\\Topo_genes_qPCR\qPCR_results.xlsx"
plant_1_pc=pd.read_excel(plant_1_pc_table, sheet_name='qPCR_primers_calibration', header=0, index_col=0)
print(plant_1_pc)


#Plot data.
def qPCR_primers_calibration(dataframe, suptitle_text, outpath):
    fig, plot_av=plt.subplots(4,3,figsize=(12,10), dpi=100)
    fig.suptitle(suptitle_text, size=15)
    
    #Prepare x axis.
    Conc_data=dataframe.loc['Concentration', :].tolist()
    
    #Plot data.
    Num_of_datasets=len(dataframe.index.tolist())-1
    print(Num_of_datasets)
    Primers_list=dataframe.index.tolist()[1:]
    print(Primers_list)
    for i in range(Num_of_datasets):
        #Points.
        Primers_pair=Primers_list[i]
        print(Primers_pair)
        qC_data=dataframe.loc[Primers_pair, :]
        print(Primers_pair, qC_data)
        print(i, int(i/3), i%3)
        plot_av[int(i/3), i%3].scatter(Conc_data, qC_data.tolist(), s=2, color='k', edgecolors='black', linewidth=0.2, alpha=1, zorder=1000) 
        
        #Linear fitting of linear data.
        idx_no_nan=np.isfinite(Conc_data) & np.isfinite(qC_data.tolist())
        print(idx_no_nan)
        fit=np.polyfit(np.array(Conc_data)[idx_no_nan], np.array(qC_data.tolist())[idx_no_nan], 1)  
        
        #fit=np.polyfit(Conc_data, qC_data.tolist(), 1)
        print(fit)
        fit_fn=np.poly1d(fit) 
        plot_av[int(i/3), i%3].plot(Conc_data, fit_fn(Conc_data), '--b', linewidth=0.5, label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3))) 
        primer_effectiveness=10**(-1/fit[0])
        plot_av[int(i/3), i%3].annotate(f"$\lambda$={round(primer_effectiveness, 2)}", xy=(0.7, 0.5), xycoords='axes fraction', size=16)
        plot_av[int(i/3), i%3].set_ylabel('Ct', size=20)
        plot_av[int(i/3), i%3].set_xlabel('Relative template concentration', size=12)
        plot_av[int(i/3), i%3].set_xticks(np.unique(Conc_data), minor=False)
        plot_av[int(i/3), i%3].set_xticklabels([1, 10, 100, 1000], rotation=0, size=12)
        #plot_av[int(i/3), i%3].set_xscale('log')
        plot_av[int(i/3), i%3].set_title(f"Primers pair {Primers_pair}", size=12)
        plot_av[int(i/3), i%3].legend(frameon=False)
        

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    plt.savefig(f'{outpath}.png', dpi=300, size=(12,13))
    plt.savefig(f'{outpath}.svg', dpi=300, size=(12,13))

    return

qPCR_primers_calibration(plant_1_pc, 'Plant qPCR primers calibration', "C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\\Topo_genes_qPCR\Plant_primers_callibration")



#################
### Arabidopsis gyrase genes expression analysis in response to Auxin treatment.
#################

#Path to qPCR data table.
data_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Topo_genes_qPCR\\qPCR_results.xlsx"
data_tab=pd.read_excel(data_table, sheet_name='qPCR_Ef_1a_norm', header=0, index_col=0)
print(data_tab)

#Compare sets of values using t-test.
def compare_sets_t_test(set_of_sets, set_names, gene_name):
    
    for i in range(len(set_of_sets)):
        
        for j in range(len(set_of_sets)):
            
            if j>i:
            
                set1=set_of_sets[i]
                set_name1=set_names[i]
                set2=set_of_sets[j]
                set_name2=set_names[j]
                
                t_test_stat=stats.ttest_ind(set1, set2)
                
                print(f'T-test for gene {gene_name} condition {set_name1} vs condition {set_name2} statistic: {t_test_stat[0]}; p-value: {t_test_stat[1]}')
            
    return


#Plot data.
def qPCR_expression_auxin(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(5,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{gyrA}$', '$\it{gyrB1}$', '$\it{gyrB2}$', '$\it{gyrB3}$']
    
    #print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8]
    
    X_coords_main=[1.9,4.9,7.9,10.9]
    
    #print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['Control mean', 'Aux 0.5 mean', 'Aux 1 mean', 'Aux 2 mean']
    Mean_Ct=dataframe.loc['GyrA', Data_points_ar].tolist() + dataframe.loc['GyrB1', Data_points_ar].tolist() + dataframe.loc['GyrB2', Data_points_ar].tolist() + dataframe.loc['GyrB3', Data_points_ar].tolist() 
    
    #print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['Control std', 'Aux 0.5 std', 'Aux 1 std', 'Aux 2 std']
    Errors_Ct=dataframe.loc['GyrA', Data_errors_ar].tolist() + dataframe.loc['GyrB1', Data_errors_ar].tolist() + dataframe.loc['GyrB2', Data_errors_ar].tolist() + dataframe.loc['GyrB3', Data_errors_ar].tolist()
    
    #print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    
    #print(len(Colors))
    
    #Prepare data for real data points.
    control_labs=['Control 1', 'Control 2', 'Control 3']
    set1_labs=['Aux 0.5 1', 'Aux 0.5 2', 'Aux 0.5 3']
    set2_labs=['Aux 1 1', 'Aux 1 2', 'Aux 1 3']
    set3_labs=['Aux 2 1', 'Aux 2 2', 'Aux 2 3']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['GyrA', control_labs].tolist(), dataframe.loc['GyrA', set1_labs].tolist(), dataframe.loc['GyrA', set2_labs].tolist(), dataframe.loc['GyrA', set3_labs].tolist(),
                 dataframe.loc['GyrB1', control_labs].tolist(), dataframe.loc['GyrB1', set1_labs].tolist(), dataframe.loc['GyrB1', set2_labs].tolist(), dataframe.loc['GyrB1', set3_labs].tolist(),
                 dataframe.loc['GyrB2', control_labs].tolist(), dataframe.loc['GyrB2', set1_labs].tolist(), dataframe.loc['GyrB2', set2_labs].tolist(), dataframe.loc['GyrB2', set3_labs].tolist(),
                 dataframe.loc['GyrB3', control_labs].tolist(), dataframe.loc['GyrB3', set1_labs].tolist(), dataframe.loc['GyrB3', set2_labs].tolist(), dataframe.loc['GyrB3', set3_labs].tolist()]
    
    #print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    #print(len(X_coords_for_points))
    #print(len(Data_points_ar))
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    #plot_av.set_yticklabels([0,0.01,0.02,0.03,0.04,0.05], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 0.085])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('Control', 'Aux 0.5', 'Aux 1', 'Aux 2'), fontsize=13, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(f'{outpath}.png', dpi=300, size=(5,3))
    plt.savefig(f'{outpath}.svg', dpi=300, size=(5,3))
    
    # Perform t-test between conditions.
    Set_names=['Control', 'Aux 0.5', 'Aux 1', 'Aux 2']
    Gene_name='GyrA'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)
    
    Gene_name='GyrB1'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)  
    
    Gene_name='GyrB2'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name) 
    
    Gene_name='GyrB3'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)     
    
    return


#Plot data.
def qPCR_expression_cytokinin(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(5,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{gyrA}$', '$\it{gyrB1}$', '$\it{gyrB2}$', '$\it{gyrB3}$']
    
    #print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8]
    
    X_coords_main=[1.9,4.9,7.9,10.9]
    
    #print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['Control mean', 'Ctk 0.25 mean', 'Ctk 0.5 mean', 'Ctk 1 mean']
    Mean_Ct=dataframe.loc['GyrA', Data_points_ar].tolist() + dataframe.loc['GyrB1', Data_points_ar].tolist() + dataframe.loc['GyrB2', Data_points_ar].tolist() + dataframe.loc['GyrB3', Data_points_ar].tolist() 
    
    #print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['Control std', 'Ctk 0.25 std', 'Ctk 0.5 std', 'Ctk 1 std']
    Errors_Ct=dataframe.loc['GyrA', Data_errors_ar].tolist() + dataframe.loc['GyrB1', Data_errors_ar].tolist() + dataframe.loc['GyrB2', Data_errors_ar].tolist() + dataframe.loc['GyrB3', Data_errors_ar].tolist()
    
    #print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    
    #print(len(Colors))
    
    #Prepare data for real data points.
    control_labs=['Control 1', 'Control 2', 'Control 3']
    set1_labs=['Ctk 0.25 1', 'Ctk 0.25 2', 'Ctk 0.25 3']
    set2_labs=['Ctk 0.5 1', 'Ctk 0.5 2', 'Ctk 0.5 3']
    set3_labs=['Ctk 1 1', 'Ctk 1 2', 'Ctk 1 3']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['GyrA', control_labs].tolist(), dataframe.loc['GyrA', set1_labs].tolist(), dataframe.loc['GyrA', set2_labs].tolist(), dataframe.loc['GyrA', set3_labs].tolist(),
                 dataframe.loc['GyrB1', control_labs].tolist(), dataframe.loc['GyrB1', set1_labs].tolist(), dataframe.loc['GyrB1', set2_labs].tolist(), dataframe.loc['GyrB1', set3_labs].tolist(),
                 dataframe.loc['GyrB2', control_labs].tolist(), dataframe.loc['GyrB2', set1_labs].tolist(), dataframe.loc['GyrB2', set2_labs].tolist(), dataframe.loc['GyrB2', set3_labs].tolist(),
                 dataframe.loc['GyrB3', control_labs].tolist(), dataframe.loc['GyrB3', set1_labs].tolist(), dataframe.loc['GyrB3', set2_labs].tolist(), dataframe.loc['GyrB3', set3_labs].tolist()]
    
    #print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    #print(len(X_coords_for_points))
    #print(len(Data_points_ar))
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    #plot_av.set_yticklabels([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 0.18])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('Control', 'Ctk 0.25', 'Ctk 0.5', 'Ctk 1'), fontsize=12, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(f'{outpath}.png', dpi=300, size=(5,3))
    plt.savefig(f'{outpath}.svg', dpi=300, size=(5,3))
    
    # Perform t-test between conditions.
    Set_names=['Control', 'Ctk 0.25', 'Ctk 0.5', 'Ctk 1']
    Gene_name='GyrA'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)
    
    Gene_name='GyrB1'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)  
    
    Gene_name='GyrB2'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name) 
    
    Gene_name='GyrB3'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)         

    return


#Plot data.
def qPCR_expression_gib(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(5,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{gyrA}$', '$\it{gyrB1}$', '$\it{gyrB2}$', '$\it{gyrB3}$']
    
    #print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8]
    
    X_coords_main=[1.9,4.9,7.9,10.9]
    
    #print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['Control mean', 'Gib 0.125 mean', 'Gib 0.25 mean', 'Gib 0.5 mean']
    Mean_Ct=dataframe.loc['GyrA', Data_points_ar].tolist() + dataframe.loc['GyrB1', Data_points_ar].tolist() + dataframe.loc['GyrB2', Data_points_ar].tolist() + dataframe.loc['GyrB3', Data_points_ar].tolist() 
    
    #print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['Control std', 'Gib 0.125 std', 'Gib 0.25 std', 'Gib 0.5 std']
    Errors_Ct=dataframe.loc['GyrA', Data_errors_ar].tolist() + dataframe.loc['GyrB1', Data_errors_ar].tolist() + dataframe.loc['GyrB2', Data_errors_ar].tolist() + dataframe.loc['GyrB3', Data_errors_ar].tolist()
    
    #print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    
    #print(len(Colors))
    
    #Prepare data for real data points.
    control_labs=['Control 1', 'Control 2', 'Control 3']
    set1_labs=['Gib 0.125 1', 'Gib 0.125 2', 'Gib 0.125 3']
    set2_labs=['Gib 0.25 1', 'Gib 0.25 2', 'Gib 0.25 3']
    set3_labs=['Gib 0.5 1', 'Gib 0.5 2', 'Gib 0.5 3']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['GyrA', control_labs].tolist(), dataframe.loc['GyrA', set1_labs].tolist(), dataframe.loc['GyrA', set2_labs].tolist(), dataframe.loc['GyrA', set3_labs].tolist(),
                 dataframe.loc['GyrB1', control_labs].tolist(), dataframe.loc['GyrB1', set1_labs].tolist(), dataframe.loc['GyrB1', set2_labs].tolist(), dataframe.loc['GyrB1', set3_labs].tolist(),
                 dataframe.loc['GyrB2', control_labs].tolist(), dataframe.loc['GyrB2', set1_labs].tolist(), dataframe.loc['GyrB2', set2_labs].tolist(), dataframe.loc['GyrB2', set3_labs].tolist(),
                 dataframe.loc['GyrB3', control_labs].tolist(), dataframe.loc['GyrB3', set1_labs].tolist(), dataframe.loc['GyrB3', set2_labs].tolist(), dataframe.loc['GyrB3', set3_labs].tolist()]
    
    #print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    #print(len(X_coords_for_points))
    #print(len(Data_points_ar))
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    #plot_av.set_yticklabels([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 0.60])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('Control', 'Gib 0.125', 'Gib 0.25', 'Gib 0.5'), fontsize=12, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(f'{outpath}.png', dpi=300, size=(5,3))
    plt.savefig(f'{outpath}.svg', dpi=300, size=(5,3))
    
    # Perform t-test between conditions.
    Set_names=['Control', 'Gib 0.125', 'Gib 0.25', 'Gib 0.5']
    Gene_name='GyrA'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)
    
    Gene_name='GyrB1'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)  
    
    Gene_name='GyrB2'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name) 
    
    Gene_name='GyrB3'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)             
    
    return


#Plot data.
def qPCR_expression_SA(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(5,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{gyrA}$', '$\it{gyrB1}$', '$\it{gyrB2}$', '$\it{gyrB3}$']
    
    #print(len(Conditions))
    
    X_coords=[1,1.6,2.2,2.8,
              4,4.6,5.2,5.8,
              7,7.6,8.2,8.8,
              10,10.6,11.2,11.8]
    
    X_coords_main=[1.9,4.9,7.9,10.9]
    
    #print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['Control mean', 'SA 0.125 mean', 'SA 0.25 mean', 'SA 0.5 mean']
    Mean_Ct=dataframe.loc['GyrA', Data_points_ar].tolist() + dataframe.loc['GyrB1', Data_points_ar].tolist() + dataframe.loc['GyrB2', Data_points_ar].tolist() + dataframe.loc['GyrB3', Data_points_ar].tolist() 
    
    #print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['Control std', 'SA 0.125 std', 'SA 0.25 std', 'SA 0.5 std']
    Errors_Ct=dataframe.loc['GyrA', Data_errors_ar].tolist() + dataframe.loc['GyrB1', Data_errors_ar].tolist() + dataframe.loc['GyrB2', Data_errors_ar].tolist() + dataframe.loc['GyrB3', Data_errors_ar].tolist()
    
    #print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',
            '#b2e69a', '#f598b8', '#f5ab87', '#89d8fa',]
    
    #print(len(Colors))
    
    #Prepare data for real data points.
    control_labs=['Control 1', 'Control 2', 'Control 3']
    set1_labs=['SA 0.125 1', 'SA 0.125 2', 'SA 0.125 3']
    set2_labs=['SA 0.25 1', 'SA 0.25 2', 'SA 0.25 3']
    set3_labs=['SA 0.5 1', 'SA 0.5 2', 'SA 0.5 3']
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['GyrA', control_labs].tolist(), dataframe.loc['GyrA', set1_labs].tolist(), dataframe.loc['GyrA', set2_labs].tolist(), dataframe.loc['GyrA', set3_labs].tolist(),
                 dataframe.loc['GyrB1', control_labs].tolist(), dataframe.loc['GyrB1', set1_labs].tolist(), dataframe.loc['GyrB1', set2_labs].tolist(), dataframe.loc['GyrB1', set3_labs].tolist(),
                 dataframe.loc['GyrB2', control_labs].tolist(), dataframe.loc['GyrB2', set1_labs].tolist(), dataframe.loc['GyrB2', set2_labs].tolist(), dataframe.loc['GyrB2', set3_labs].tolist(),
                 dataframe.loc['GyrB3', control_labs].tolist(), dataframe.loc['GyrB3', set1_labs].tolist(), dataframe.loc['GyrB3', set2_labs].tolist(), dataframe.loc['GyrB3', set3_labs].tolist()]
    
    #print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    #print(len(X_coords_for_points))
    #print(len(Data_points_ar))
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    #plot_av.set_yticklabels([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 0.2])
    
    plt.legend((Bars[0],Bars[1],Bars[2],Bars[3]), ('Control', 'SA 0.125', 'SA 0.25', 'SA 0.5'), fontsize=12, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(f'{outpath}.png', dpi=300, size=(5,3))
    plt.savefig(f'{outpath}.svg', dpi=300, size=(5,3))
    
    # Perform t-test between conditions.
    Set_names=['Control', 'SA 0.125', 'SA 0.25', 'SA 0.5']
    Gene_name='GyrA'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)
    
    Gene_name='GyrB1'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)  
    
    Gene_name='GyrB2'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name) 
    
    Gene_name='GyrB3'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist(), dataframe.loc[Gene_name, set2_labs].tolist(), dataframe.loc[Gene_name, set3_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)    
    
    return


#Plot data.
def qPCR_expression_dark(dataframe, outpath):
    
    ###
    ##Plot all FEs.
    ###
    
    fig, plot_av=plt.subplots(1,1,figsize=(4,3), dpi=100)
    
    #Prepare x axis.
    Conditions=['$\it{gyrA}$', '$\it{gyrB1}$', '$\it{gyrB2}$', '$\it{gyrB3}$']
    
    #print(len(Conditions))
    
    X_coords=[1,1.6,
              2.8,3.4,
              4.6,5.2,
              6.4,7.0]
    
    X_coords_main=[1.3,3.1,4.9,6.7]
    
    #print(len(X_coords))
    
    #Prepare data for bars.
    Data_points_ar=['Control mean', 'Dark mean']
    Mean_Ct=dataframe.loc['GyrA', Data_points_ar].tolist() + dataframe.loc['GyrB1', Data_points_ar].tolist() + dataframe.loc['GyrB2', Data_points_ar].tolist() + dataframe.loc['GyrB3', Data_points_ar].tolist() 
    
    #print(len(Mean_Ct))
    
    #Prepare data for error bars (precomputed standard deviation).
    Data_errors_ar=['Control std', 'Dark std']
    Errors_Ct=dataframe.loc['GyrA', Data_errors_ar].tolist() + dataframe.loc['GyrB1', Data_errors_ar].tolist() + dataframe.loc['GyrB2', Data_errors_ar].tolist() + dataframe.loc['GyrB3', Data_errors_ar].tolist()
    
    #print(len(Errors_Ct))
    
    #Set colors for bars.
    Colors=['#b2e69a', '#f598b8',
            '#b2e69a', '#f598b8',
            '#b2e69a', '#f598b8',
            '#b2e69a', '#f598b8',]
    
    #print(len(Colors))
    
    #Prepare data for real data points.
    control_labs=['Control 1', 'Control 2', 'Control 3']
    set1_labs=['Dark 1', 'Dark 2', 'Dark 3']  
    
    #For some reason plotting of a list of lists doesn't work if sublists are of unequal size. E.g. [x1, x2], [[y1,y2],[y3,y4]] works fine and [x1, x2], [[y1],[y3,y4]] does not.
    Data_points=[dataframe.loc['GyrA', control_labs].tolist(), dataframe.loc['GyrA', set1_labs].tolist(),
                 dataframe.loc['GyrB1', control_labs].tolist(), dataframe.loc['GyrB1', set1_labs].tolist(), 
                 dataframe.loc['GyrB2', control_labs].tolist(), dataframe.loc['GyrB2', set1_labs].tolist(), 
                 dataframe.loc['GyrB3', control_labs].tolist(), dataframe.loc['GyrB3', set1_labs].tolist(),]
    
    #print(len(Data_points))
    
    X_coords_for_points=[]
    Data_points_ar=[]
    
    for i in range(len(Data_points)):
        X_coords_for_points=X_coords_for_points+[X_coords[i]]*len(Data_points[i])
        Data_points_ar=Data_points_ar+Data_points[i]
    
    #print(len(X_coords_for_points))
    #print(len(Data_points_ar))
    
    #Plot data.
    Bars=plot_av.bar(X_coords, Mean_Ct, yerr=Errors_Ct, error_kw=dict(lw=1, capsize=3, capthick=1), align='center', width=0.6, color=Colors, edgecolor='k', linewidth=0.6)
    plot_av.plot(X_coords_for_points, Data_points_ar, 'ko', markersize=1) 
    plot_av.set_ylabel('Normalized value', size=17)
    plot_av.set_xticks(X_coords_main)
    plot_av.set_xticklabels(Conditions, rotation=0, size=14)  
    #plot_av.set_yticklabels([0,0.01,0.02,0.03,0.04,0.05], rotation=0, size=12)
    plot_av.tick_params(axis='x', which='major', pad=5)
    plot_av.set_ylim([0, 0.075])
    
    plt.legend((Bars[0],Bars[1]), ('Control', 'Dark'), fontsize=14, ncol=4, frameon=False, markerscale=2, handlelength=0.7, handletextpad=0.3, columnspacing=0.7)
    plt.tight_layout()
    
    plt.show()    
    plt.savefig(f'{outpath}.png', dpi=300, size=(4,3))
    plt.savefig(f'{outpath}.svg', dpi=300, size=(4,3))
    
    # Perform t-test between conditions.
    Set_names=['Control', 'Dark']
    Gene_name='GyrA'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)
    
    Gene_name='GyrB1'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)  
    
    Gene_name='GyrB2'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name) 
    
    Gene_name='GyrB3'
    Set_of_sets=[dataframe.loc[Gene_name, control_labs].tolist(), dataframe.loc[Gene_name, set1_labs].tolist()]
    compare_sets_t_test(Set_of_sets, Set_names, Gene_name)     
    
    return

qPCR_expression_auxin(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Topo_genes_qPCR\Test\Arabidopsis_gyrase_expression_Aux_EF_1a")
qPCR_expression_cytokinin(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Topo_genes_qPCR\Test\Arabidopsis_gyrase_expression_Ctk_EF_1a")
qPCR_expression_gib(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Topo_genes_qPCR\Test\Arabidopsis_gyrase_expression_Gib_EF_1a")
qPCR_expression_SA(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Topo_genes_qPCR\Test\Arabidopsis_gyrase_expression_SA_EF_1a")
qPCR_expression_dark(data_tab, "C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Topo_genes_qPCR\Test\Arabidopsis_gyrase_expression_Dark_EF_1a")