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



def nuc_div(pop,pos):
    """
    Estimate Nucleotide diversity.

    
    Input: population genotype array data for each population (calculated separately )
    pop = [[[0,0],[0,0]],[[1,1],[1,1]]] >> represents 2 variants and 2 samples in each snp
    pos = [1,2,3] >> represents location of the snp
    converted to scikit allel AllelCount array data structure

    Returns:	
        pi, Nucleotide diversity. 

    """
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    pi= allel.sequence_diversity(pos=pos, ac=ac)
    
    return round(pi,3)


def haplotype_div(pop):
    """
    Estimate nucleotide diversity within a given region, which is the average proportion of sites
    that differ between randomly chosen pairs of chromosomes.


    Input: population genotype array data for each population(calculated separately )
    pop = [[[0,0],[0,0]],[[1,1],[1,1]]] >> represents 2 variants and 2 samples in each snp
    converted to scikit allel Haplotype array data structure

    Returns:
        hd,Haplotype diversity.
    
    """
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    h = pop.to_haplotypes()
    hd= allel.haplotype_diversity(h=h)
    
    return round(hd,3)


def tajima_d(pop, pos):
    """
    Calculate the value of Tajima's D over a given region.


    Input: population genotype array data for each population(calculated separately )
    pop = [[[0,0],[0,0]],[[1,1],[1,1]]] >> represents 2 variants and 2 samples in each snp
    pos = [1,2] >> represents genome location of the 2 variants above
    converted to scikit allel AllelCount array data structure

    Returns:
        D, Tajima's D
    
    """
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    D= allel.tajima_d(ac=ac, pos=pos )
    
    return round(D,3)



def hudson_fst(pop1,pop2):
    """
    Need two population genotype array
    pop1 : [[[0,0],[0,0]],[[1,1],[1,1]]] >> represents 2 variants and 2 samples in each snp from the first population
    pop2 : [[[1,1],[0,1]],[[0,0],[1,1]]] >> represents 2 variants and 2 samples in each snp from the other population

    converts the above data into a scikit AllelCount array data structure, used to calculate fst
    num and den : output numerator and denominator arrays for each variant, which are averaged and returned
    """
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
    e.g. {'GBR': [[[0,0],[0,1]],[[0,0],[1,1]], 'JPT':[[[0,0],[0,1]],[[0,0],[1,1]]} >> Representing 2 snps and 2 samples from JPT and GBR population

    performs permutations to create all pairwise comparisions in a tuple
    
    Input:
    [GBR,JPT,YRI] 
    
    Outputs:
    [('GBR','JPT'),('GBR','YRI'),('JPT','GBR'),('JPT','YRI'),('YRI','GBR'), ('YRI', 'JPT')]. 
    
    For each permutation, it calculates the fst by using the items in tuple to only call the 
    genotype data dictionary that is specific to the population.  

    Returns:
    A list with first element indicating the pairs compared and the second indication the fst value of those comparisions
    i.e. [('GBR', 'JPT'), 0.3] >> represents the fst value of GBR
    """
    combos = itertools.permutations(paris,2)
    lst = []
    for i in combos:
        pair = f'{i[0]}:{i[1]}'
        # passes the two genotype data arrays each represented by their keys 
        fst = hudson_fst(gt_dict[i[0]],gt_dict[i[1]])
        lst.append([pair,fst])
    
    return lst








def overall_stats_gtd(all_pops_array):
    """ 
    For calculating the stats using scikit allel using genotype data and combining 
    all the individuals into a single super population, by using numpy hstack and iterating
    over the elements in all_pop_array, each representing the population selected

    """
    for i,ele in enumerate(all_pops_array):
        ele = np.array(ele)
        if i == 0:
            # creates and assigns the array for the first element
            all_pop_arr = ele
        else:
            # stacks the arrays after the first element whilst maintaing the dimensions
            all_pop_arr = np.hstack((all_pop_arr,ele))
    return all_pop_arr

    

def overall_stats_cts(all_pop):
    """ 
    For calculating observed Homozygosity using count data and combining 
    all the individuals into a single super population
    """
    full_pop_data = []
    if len(full_pop_data) == 0:
        frq_lst = all_pop[0]
        # Update the counts from a str to int type for calulations
        for dt in frq_lst:
            for k,v in dt.items():
                dt.update({k: int(v)})
        for i in frq_lst:
            full_pop_data.append(i)
    frq_lst = all_pop[1:]
    # Add all the counts to their respective key from all population, for each snp.
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
    """
    Estimate nucleotide diversity in windows over range of snps, selected by the user.
    
    First element of the snp location as start and, the last position as the stop to create a range
    pop is the genotype array data for each population i.e. [[[0,0],[0,0]],[[1,1],[1,1]]] >>  represents 2 variants and 2 samples in each snp
    bin_size refers to the window size which is user selected
    step_size refers to the overlapping windows also specified by the user

    Returns:
        Pi values in an array, each value corresponding to the window  
    """   
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    win_pi, windows, n_bases, counts = allel.windowed_diversity(pos,ac,size=bin_size,step=step_size)
    return win_pi



def win_tajima_d(positions,pop,bin_size=100,step_size=None):
    """
    Calculates Tajimas D in windows over range of snps, selected by the user.

    First element of the snp location as start and, the last position as the stop to create a range
    pop is the genotype array data for each population i.e. [[[0,0],[0,0]],[[1,1],[1,1]]] >>  represents 2 variants and 2 samples in each snp
    bin_size refers to the window size which is user selected
    step_size refers to the overlapping windows also specified by the user

    Returns:
        D values in an array, each value corresponding to the window
    """
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    ac = pop.count_alleles()
    win_tajima_D, windows, counts = allel.windowed_tajima_d(pos=pos,ac=ac,size=bin_size,step=step_size,min_sites=2)
    return win_tajima_D



def win_haplotype_div(positions,pop,bin_size=100,step_size=None):
    """
    Estimate Haplotype diversity in windows over range of snps, selected by the user.

    First element of the snp location as start and, the last position as the stop to create a range
    pop is the genotype array data for each population i.e. [[[0,0],[0,0]],[[1,1],[1,1]]] >>  represents 2 variants and 2 samples in each snp
    bin_size refers to the window size which is user selected
    step_size refers to the overlapping windows also specified by the user

    Returns:
        hd values in an array, each value corresponding to the window
    """
    pos = np.array(positions)
    pop_gt = np.array(pop)
    pop = allel.GenotypeArray(pop_gt)
    pop_hap = pop.to_haplotypes()
    windowed_hap, windows, counts = allel.windowed_statistic(pos=pos,values=pop_hap,statistic=allel.haplotype_diversity,size=bin_size, step=step_size)
    return windowed_hap



def win_hudson_fst(positions,pop1,pop2,bin_size=100,step_size=None):
    """
    Calculate FST in windows over range of snps, selected by the user.

    'USES get_win_fstat() to get the combinations'

    First element of the snp location as start and, the last position as the stop to create a range
    pop1 is the genotype array data for one population i.e. [[[0,0],[0,0]],[[1,1],[1,1]]] >>  represents 2 variants and 2 samples in each snp
    pop2 is the genotype array data for the other population i.e. [[[0,0],[0,0]],[[1,1],[1,1]]] >>  represents 2 variants and 2 samples in each snp
    bin_size refers to the window size which is user selected
    step_size refers to the overlapping windows also specified by the user

    Returns:
        fst values in an array, each value corresponding to the window for each pair
    """
    pop1_gt = np.array(pop1)
    pop2_gt = np.array(pop2)
    pop1 = allel.GenotypeArray(pop1_gt)
    pop2 = allel.GenotypeArray(pop2_gt)
    ac1 = pop1.count_alleles()
    ac2 = pop2.count_alleles()

    win_fst, windows, counts = allel.windowed_hudson_fst(positions,ac1, ac2,size=bin_size,step=step_size)

    return win_fst



def get_win_fstat(paris,gt_dict,pos,bin_size,step_size):
    """ Windowed hudson Fst using permutations to get pair-wise combinations of populations similar to 
        the get_fstat() and uses the  win_hudson_fst() to calculate windowed fst for each combination
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
    """ Creates the x axis  by gett the average of the window size"""
    windows = allel.position_windows(pos=pos, size=size, step=step, start=None, stop=None)
    x = np.asarray(windows).mean(axis=1)
    return x


def plot_win_taj_d(TD, position, num_pops):

    # TD consists of dictionary with keys represent the population and the windowed Tajima D values
    # graph containing all populations stored as a dictionary with keys representing the population
    # positions are the values calculated from avg_win()
        
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
    # HD consists of dictionary with keys represent the population and the windowed Haplotype diversity values
    # e.g. {'GBR': [0.1,0.2], 'JPT':[0.3,0.4]} >> represents GBR and JPT population haplotype diversity in 2 windows 
    # graph containing all populations stored as a dictionary with keys representing the population
    # positions are the values calculated from avg_win()
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

    # ND consists of dictionary with keys represent the population and the windowed Nucleotide diversity values
    # e.g. {'GBR': [0.1,0.2], 'JPT':[0.3,0.4]} >> represents GBR and JPT population nucleotide diversity in 2 windows
    # graph containing all populations stored as a dictionary with keys representing the population
    # positions are the values calculated from avg_win()
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

    # pop_FST consists of dictionary with keys represent each pariwise comparisons and the FST Nucleotide diversity values
    # e.g. {'GBR-JPT': [0.1,0.2], 'GBR-MXL':[0.3,0.4]} >> represents GBR-JPT pairs FST values and 'GBR-MXL' population nucleotide diversity in 2 windows
    # Plots are split based on populations, i.e. all GBR in one plot and all JPT comparisions are in another plot
    # graph containing all populations stored as a dictionary with keys representing the population
    # positions are the values calculated from avg_win()
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
