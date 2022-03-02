from tracemalloc import start, stop
import plotly.graph_objects as go
import plotly
import json
import plotly.express as px
import allel
import ast
import pandas as pd
import numpy as np
import itertools


def Homozygosity(freq_data):
    ''' Estimates a score for homozygosity using direct count method as explained in Sabatti & Risch, 2002. 
    Input a list of dictionaries retrieved from search query. 
    Each item should correspond to a single rs ID.

    Parameters:
        freq_data (dict): genotype count data. Keys correspond to hom_ref, het and hom_alt
    Return: 
        Observed homozygosity (rounded to 3 d.p)
    '''
    total_hom = 0
    total_gen = 0
    for x in freq_data:
        total_hom += int(x['hom_ref']) + int(x['hom_alt'])
        total_gen += int(x['hom_ref']) + int(x['hom_alt']) + int(x['het'])
    obs_homozygosity = total_hom/total_gen
    return round(obs_homozygosity,3)


# Using scikit-allel - windows determined by number of variants. will need to be plotted with
# rs IDs on x-axis to work correctly. Exception being nucleotide diversity.
# bin_size and step_size are user submitted, and pop refers to a single population chosen. functions
# would need to be run for each population if multiple populations are selected. 
def nuc_div(pop,pos):
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    pi= allel.sequence_diversity(pos=pos, ac=ac)
    
    return round(pi,3)


def haplotype_div(pop):
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    h = pop.to_haplotypes()
    hd= allel.haplotype_diversity(h=h)
    
    return round(hd,3)


def tajima_d(pop, pos):
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    D= allel.tajima_d(ac=ac, pos=pos )
    
    return round(D,3)


def hudson_fst(pop1,pop2):
    pop1_gt = np.array(pop1)
    pop2_gt = np.array(pop2)
    pop1 = allel.GenotypeArray(pop1_gt)
    pop2 = allel.GenotypeArray(pop2_gt)
    ac1 = pop1.count_alleles()
    ac2 = pop2.count_alleles()
    num, den = allel.hudson_fst(ac1, ac2)
    fst = np.sum(num) / np.sum(den)
    return round(fst,3)


def get_main_stats(pop,freq_data,pos,stats):
    """ Bundles other basic statistics into a single functiona and also 
    checks for the user input before calculation.

    Parameters:
        pop (array): genotype array data
        freq_data (dict): genotype count data. Keys correspond to hom_ref, het and hom_alt
        pos (list): Positions for SNPs
        stats (list): Stats to be calculated
    Return:
        homo (float):
        nd ():
        hap_div ():
        taj_d ():
    
    """
    if 'Homozygosity' in stats:
        homo = Homozygosity(freq_data)
    else:
        homo = 'Not calculated'

    if 'Nucleotide Diversity' in stats:
        nd = nuc_div(pop,pos=pos)
    else:
        nd = 'Not calculated'
    
    if 'Haplotype Diversity' in stats:
        hap_div = haplotype_div(pop)
    else:
        hap_div = 'Not calculated'

    if 'Tajimas D' in stats:
        taj_d = tajima_d(pop,pos=pos)
    else:
        taj_d = 'Not calculated'

    return homo,nd,hap_div,taj_d


def get_fstat(paris,gt_dict):
    """ First get the genotype data in a dictionary format with their populations code as the key.
    performs permutations to create all pairwise comparisions in a tuple e.g.('GBR','JPT'). 
    For each permutation, it calculates the fst by using the items in tuple to only call the 
    genotype data that is specific to the population.  
    """
    combos = itertools.permutations(paris,2)
    lst = []
    for i in combos:
        pair = f'{i[0]}:{i[1]}'
        fst = hudson_fst(gt_dict[i[0]],gt_dict[i[1]])
        lst.append([pair,fst])
    
    return lst








def overall_stats_gtd(all_pops_array):
    """ For calculating the stats using scikit allel using genotype data and combining 
    all the individuals into a single super population 
    """
    for i,ele in enumerate(all_pops_array):
        ele = np.array(ele)
        if i == 0:
            all_pop_arr = ele
        else:
            all_pop_arr = np.hstack((all_pop_arr,ele))
    return all_pop_arr

    

def overall_stats_cts(all_pop):
    """ For calculating observed Homozygosity using count data and combining 
    all the individuals into a single super population
    """
    full_pop_data = []
    if len(full_pop_data) == 0:
        frq_lst = all_pop[0]
        # In place update of the dictionary from
        for dt in frq_lst:
            for k,v in dt.items():
                dt.update({k: int(v)})
        for i in frq_lst:
            full_pop_data.append(i)
    frq_lst = all_pop[1:]
    for i in frq_lst:
        for main,sub in zip(full_pop_data,i):
            main['hom_ref'] += int(sub['hom_ref'])
            main['het'] += int(sub['het'])
            main['hom_alt'] += int(sub['hom_alt'])

    return full_pop_data

  

############################################################################################################################################
######################################################    WINDOWED STATS ###################################################################
############################################################################################################################################
    

def win_nuc_div(positions,pop,bin_size=100,step_size=None):
    # input list of queries retrieved from the results page, window size and step size. 
    # window size and step size refers to number of nucleotides. produces a list of 4 arrays
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    win_pi, windows, n_bases, counts = allel.windowed_diversity(pos,ac,size=bin_size,step=step_size)
    return win_pi


def win_tajima_d(positions,pop,bin_size=100,step_size=None):
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    win_tajima_D, windows, counts = allel.windowed_tajima_d(pos=pos,ac=ac,size=bin_size,step=step_size,min_sites=2)
    return win_tajima_D


def win_haplotype_div(positions,pop,bin_size=100,step_size=None):
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    pop_hap = pop.to_haplotypes()
    windowed_hap, windows, counts = allel.windowed_statistic(pos=pos,values=pop_hap,statistic=allel.haplotype_diversity,size=bin_size, step=step_size)
    return windowed_hap


def win_hudson_fst(positions,pop1,pop2,bin_size=100,step_size=None):
    pop1_gt = np.array(pop1)
    pop2_gt = np.array(pop2)
    pop1 = allel.GenotypeArray(pop1_gt)
    pop2 = allel.GenotypeArray(pop2_gt)
    ac1 = pop1.count_alleles()
    ac2 = pop2.count_alleles()

    win_fst, windows, counts = allel.windowed_hudson_fst(positions,ac1, ac2,size=bin_size,step=step_size)

    return win_fst



def get_win_fstat(paris,gt_dict,pos,bin_size,step_size):
    """ Windowed hudson Fst using permutations to get pair-wise combinations of populations
    """
    combos = itertools.permutations(paris,2)
    lst = []
    for i in combos:
        pair = f'{i[0]}:{i[1]}'
        fst= win_hudson_fst(positions=pos,pop1 = gt_dict[i[0]], pop2= gt_dict[i[1]],bin_size=bin_size,step_size=step_size)
        lst.append([pair,fst])
    
    return lst





############################################################################################################################################
########################################################  Creating Plots ###################################################################
############################################################################################################################################


def avg_win(pos,size,step):
    windows = allel.position_windows(pos=pos, size=size, step=step, start=None, stop=None)
    x = np.asarray(windows).mean(axis=1)
    return x


def plot_win_taj_d(TD, position, num_pops):

        # graph containing all populations without buttons 
        fig_1 = go.Figure()
        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_1.add_trace(obj)

        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")


        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)

        return graph1JSON





def plot_win_hap(HD, position,num_pops):

    # graph containing all populations without buttons
    fig_2 = go.Figure()
    for x in HD.items():
        obj = go.Scatter(name=x[0], x=position, y=x[1])
        fig_2.add_trace(obj)

    fig_2.update_layout(
        xaxis_title='Position (Base pairs)',
        yaxis_title="Haploid Diversity")

 

    graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)




    return graph2JSON




def plot_nuc_div(ND, position,num_pops):

    # graph containing all populations without buttons
    fig_3 = go.Figure()
    for x in ND.items():
        obj = go.Scatter(name=x[0], x=position, y=x[1])
        fig_3.add_trace(obj)

    fig_3.update_layout(
        xaxis_title='Position (Base pairs)',
        yaxis_title="Nucleotide Diversity")

        
    graph3JSON = json.dumps(fig_3, cls=plotly.utils.PlotlyJSONEncoder)



    return graph3JSON




def plot_win_FST(pop_FST, position):
    # Get the 
    fig_7 = go.Figure()
    buttons = []
    i = 0

    # iterating through dictionary and adding each population
    for key,value in pop_FST.items():
        obj = go.Scatter(name=key, x=position, y=value)
        fig_7.add_trace(obj)

    # args is a list of booleans that tells the buttons which trace to show on click
        args = [False] * len(pop_FST)
        args[i] = True

    # create button object for each pop
        button = dict(label=key,
                        method="update",
                        args=[{"visible": args}])

    # add button to list
        buttons.append(button)
        i += 1

    # add buttons
    fig_7.update_layout(
        updatemenus=[
            dict(
                type="dropdown",
                direction="down",
                buttons=buttons)
        ])
    # add axis names
    fig_7.update_layout(
        xaxis_title='Position (Base pairs)',
        yaxis_title="FST")

    fstJSON = json.dumps(fig_7, cls=plotly.utils.PlotlyJSONEncoder)



    return fstJSON
