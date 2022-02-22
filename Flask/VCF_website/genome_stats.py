# Stats script
from scipy.stats import bartlett, chisquare
import allel
import ast
import pandas as pd
import numpy as np
import itertools
def Homozygosity(freq_data):
    # input a list of dictionaries retrieved from search query. Each item should
    # correspond to a single rs ID. 
    total_hom = 0
    total_gen = 0
    for x in freq_data:
        total_hom += int(x['hom_ref']) + int(x['hom_alt'])
        total_gen += int(x['hom_ref']) + int(x['hom_alt']) + int(x['het'])
    obs_homozygosity = total_hom/total_gen
    return obs_homozygosity




# Using scikit-allel - windows determined by number of variants. will need to be plotted with
# rs IDs on x-axis to work correctly. Exception being nucleotide diversity.
# bin_size and step_size are user submitted, and pop refers to a single population chosen. functions
# would need to be run for each population if multiple populations are selected. 


def nuc_div(pop):
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    h = pop.to_haplotypes()
    pi= allel.haplotype_diversity(h=h)
    
    return pi



def haplotype_div(pop):
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    h = pop.to_haplotypes()
    hd= allel.haplotype_diversity(h=h)
    
    return hd



def tajima_d(pop):
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    D= allel.tajima_d(ac=ac)
    
    return D





def hudson_fst(pop1,pop2):
    pop1_gt = np.array(pop1)
    pop2_gt = np.array(pop2)
    pop1 = allel.GenotypeArray(pop1_gt)
    pop2 = allel.GenotypeArray(pop2_gt)
    ac1 = pop1.count_alleles()
    ac2 = pop2.count_alleles()
    num, den = allel.hudson_fst(ac1, ac2)
    fst = np.sum(num) / np.sum(den)
    return fst

def get_fstat(paris,gt_dict):
    combos = itertools.permutations(paris,2)
    lst = []
    for i in combos:
        pair = f'{i[0]}:{i[1]}'
        fst = hudson_fst(gt_dict[i[0]],gt_dict[i[1]])
        lst.append([pair,fst])
    
    return lst





def get_main_stats(pop,freq_data):
    homo = Homozygosity(freq_data)
    nd = nuc_div(pop)
    taj_d = tajima_d(pop)
    hap_div = haplotype_div(pop)
    return homo,nd,hap_div,taj_d


def overall_stats(all_pops_array):
    for i,ele in enumerate(all_pops_array):
        ele = np.array(ele)
        if i == 0:
            all_pop_arr = ele
        else:
            all_pop_arr = np.hstack((all_pop_arr,ele))

    

    
    
    



def win_tajima_d(positions,pop,bin_size=100,step_size=None):
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    win_tajima_D, windows, counts = allel.windowed_tajima_d(pos=pos,ac=ac,size=bin_size,step=step_size)
    
    return win_tajima_D, windows, counts




def moving_haplotype_div(pop,bin_size=100,step_size=None):
    # input list of queries retrieved from the results page, window size and step size. 
    # window size refers to number of variants.
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    pop_hap = pop.to_haplotypes()
    moving_hap = allel.moving_haplotype_diversity(h=pop_hap,size=bin_size,step=step_size)
    return moving_hap





def win_nuc_div(positions,pop,bin_size=100,step_size=None):
    # input list of queries retrieved from the results page, window size and step size. 
    # window size and step size refers to number of nucleotides. produces a list of 4 arrays
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    win_pi, windows, n_bases, counts = allel.windowed_diversity(pos,ac,size=bin_size,step=step_size)
    
    return win_pi, windows, n_bases, counts


def win_hudson_fst(positions,pop1,pop2,bin_size=100,step_size=None):
    pop1_gt = np.array(pop1)
    pop2_gt = np.array(pop2)
    pop1 = allel.GenotypeArray(pop1_gt)
    pop2 = allel.GenotypeArray(pop2_gt)
    ac1 = pop1.count_alleles()
    ac2 = pop2.count_alleles()
    win_fst, win, counts = allel.windowed_hudson_fst(positions,ac1, ac2,bin_size=bin_size,step_size=step_size)
    return win_fst, win, counts











# Experimental - window size determined by chromosomal position. Window_size and steps is user 
# submitted (in bp). Sliding_window and bins functions will need to be run once. 

def sliding_window(results,window_size):
    # SNPs will need to be divided into windows based on positions, meaning that windows may not 
    # neccesarily contain the same number of SNPs
    positions = [i for i in results.pos] # grab positions for each SNP
    x = positions[0] # assigning first starting window position
    end_pos = positions[-1]
    window = [x, x + window_size] # assigning first window (start and end)
    windows = [window]
    while window[1] <= end_pos:
        x = window[1] + 1
        window = [x,x+window_size]
        windows.append(window)
    return windows

def bins(windows,results):
    windows = pd.DataFrame(windows)
    windows.columns = ['start','end']
    snp_count = []
    for row in windows.iterrows():
        count = 0
        for i in results.pos:
            if i in range [windows['start'],windows['end']]: # checking if snp is in window
                count += 1
        snp_count.append(count)
    windows['SNP'] = snp_count # dataframe include window range and number of SNPs
    bins = windows
    return bins




# unfinished
# def haplotype_div(pop,windows,step_size=None):
#     pop_array = []
#     for x in pop:
#         lst = ast.literal_eval(x)
#         pop_array.append(lst)
#     return pop_array
# # unfinished


# def tajima_d(pop,windows,step_size=None):
#     pop_array = []
#     for x in pop:
#         lst = ast.literal_eval(x)
#         pop_array.append(lst)
#     return pop_array




#############################################################################################################################################
########################################################## Lavanyas Script###################################################################
#############################################################################################################################################






#  Calculating observed and expected

def h_exp_obs_Aa(AA, Aa, aa):
    
    'Calculates expected and observed heterzygosity'
    
    p = ((2 * AA) + Aa)/ (2 * (AA + Aa + aa))
    q = ((2 * aa) + Aa)/ (2 * (AA + Aa + aa))
    
    #2pq
    twice_p_q = 2 * p * q
    
    #Calculating Observed
    Observed = (Aa *2) /(2 * (AA + Aa + aa))
 
    return Observed, twice_p_q



def H_W_Equlibrium(AA, Aa, aa):
    
    'Calculates whether there is a significant difference between the observed and expected data.'
    
    data = [AA, Aa, aa]
    chi, p= chisquare(data)
    return p


def equal_variance(O, E):     
    stat, p = bartlett(O, E)
    return p

def obs_vs_het_chi (pop_freq):
#Observed and expected heterzygosity & chi squared

    hom_ref =[] #12:14        # didnt know how to iterate using dictionary so used indices
    het = [] # 22:24
    hom_alt = [] #37:38
    rs = [] # rs values

    for x in pop_freq:
        hom_ref.append(x['hom_ref'])
        het.append(x['het'])
        hom_alt.append(x['hom_alt'])
        
    # GBR_rs = df['ID']
    # for x in GBR_rs:
    #     rs.append(x)
        
    freq_data = list(zip(hom_ref, het, hom_alt))  #getting the data
    # GBR = list(zip(rs, GBR_data))
    #print(GBR)

    # rs_value = []
    C = [] #observed and expected
    E = [] #H_W

    for x in freq_data:
        try:
            # rs_value.append(x[0])
            Calculations = h_exp_obs_Aa(AA = int(x[0]), Aa = int(x[1]), aa = int(x[2]))
            C.append(Calculations)
            Equ = H_W_Equlibrium(AA = int(x[0]), Aa = int(x[1]), aa = int(x[2]))
            E.append(Equ)
        except:
            Calculations = 'N/A' 
            Equ = 'N/A'    
            C.append(Calculations)
            E.append(Equ)
    O_E = list(zip(C, E, ))
    return O_E



