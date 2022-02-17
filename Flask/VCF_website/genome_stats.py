# Stats script
import allel
import ast
import pandas as pd
import numpy as np

def Homozygosity(population):
    # input a list of dictionaries retrieved from search query. Each item should
    # correspond to a single rs ID. 
    total_hom = 0
    total_gen = 0
    for x in population:
        geno_freq = ast.literal_eval(x['geno_freq'])
        total_hom += int(geno_freq['hom_ref']) + int(geno_freq['hom_alt'])
        total_gen += int(geno_freq['hom_ref']) + int(geno_freq['hom_alt']) + int(geno_freq['het'])
    obs_homozygosity = total_hom/total_gen
    return obs_homozygosity

# Using scikit-allel - windows determined by number of variants. will need to be plotted with
# rs IDs on x-axis to work correctly. Exception being nucleotide diversity.
# bin_size and step_size are user submitted, and pop refers to a single population chosen. functions
# would need to be run for each population if multiple populations are selected. 

def haplotype_div(pop,bin_size=100,step_size=None):
    # input list of queries retrieved from the results page, window size and step size. 
    # window size refers to number of variants.
    pop_array = []
    for x in pop:
        pop_array.append(ast.literal_eval(x['genotypes']))
    pop = allel.GenotypeArray(np.array(pop_array))
    pop_hap = pop.to_haplotypes()
    results = allel.moving_haplotype_diversity(h=pop_hap,size=bin_size,step=step_size)
    return results

def tajima_d(positions,pop,bin_size=100,step_size=None):
    pos = np.array(positions)
    pop_array = []
    for x in pop:
        pop_array.append(ast.literal_eval(x['genotypes']))
    pop = allel.GenotypeArray(pop_array)
    ac = pop.count_alleles()
    tajima_D, windows, counts = allel.windowed_tajima_d(pos=pos,ac=ac,size=bin_size,step=step_size)
    
    return tajima_D,windows,counts

def nucleotide_div(positions,pop,bin_size=100,step_size=None):
    # input list of queries retrieved from the results page, window size and step size. 
    # window size and step size refers to number of nucleotides. produces a list of 4 arrays
    pos = np.array(positions)
    pop_array = []
    for x in pop:
        pop_array.append(ast.literal_eval(x['genotypes']))
    pop = positions
    pop = allel.GenotypeArray(pop_array)
    ac = pop.count_alleles()
    pi, windows, n_bases, counts = allel.windowed_diversity(pos,ac,size=bin_size,step=step_size)
    
    return pi,windows,n_bases,counts


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