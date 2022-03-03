import string
import pandas as pd
import allel
import zarr
import numcodecs
import numpy as np
import sys
import ast


# Reading in vcf data
vcf_path = 'snp.vcf.gz'
zarr_path = 'snp.zarr'
#allel.vcf_to_zarr(vcf_path, zarr_path, group='22', fields='*', log=sys.stdout, overwrite=True)

# setting callset to zarr file
callset = zarr.open_group(zarr_path, mode='r')

# calling variants and obtaining dataframe of id, positions, reference and alternate alleles
variants = allel.VariantChunkedTable(callset['22/variants'],names=['POS','ID','REF','ALT','DP','numalt','QUAL','FILTER_PASS'],index='POS')
variants['FILTER_PASS'][:]
id = variants['ID'][:]
pos = variants['POS'][:]
ref = variants['REF'][:]
alt = variants['ALT'][:]
query_df = np.column_stack((id,pos,ref,alt))
query_df = pd.DataFrame(query_df)
query_df = query_df[[0,1,2,3]]
query_df.columns = ('ID','POS','REF','ALT')
query_df.to_csv('query_new.tsv',sep='\t',index=False)
query_df


# Assigning calldata variable
calldata = callset['22/calldata']
list(calldata.keys())

# reading in genotype data and sample information
genotypes = allel.GenotypeChunkedArray(calldata['GT'])
genotypes
samples = pd.read_table('sample_final.tsv', sep='\t')

# generating dictionary of populations and indicies
subpops = {
    'all': list(range(len(samples))),
    'MXL': samples[samples.population_code == 'MXL'].index.tolist(),
    'GBR': samples[samples.population_code == 'GBR'].index.tolist(),
    'JPT': samples[samples.population_code == 'JPT'].index.tolist(),
    'PJL': samples[samples.population_code == 'PJL'].index.tolist(),
    'YRI': samples[samples.population_code == 'YRI'].index.tolist(),
}
subpops['MXL']

# allele counts for different populations
ac_subpops = genotypes.count_alleles_subpops(subpops,max_allele=1)
MXL_ac = ac_subpops['MXL']
GBR_ac = ac_subpops['GBR']
JPT_ac = ac_subpops['JPT']
PJL_ac = ac_subpops['PJL']
YRI_ac = ac_subpops['YRI']
all_ac = ac_subpops['all']
ac_subpops


# subsetting genotype information for each population
MXL_sample = samples[samples.population_code == 'MXL'].index
MXL_g = genotypes.subset(sel0=None, sel1=MXL_sample)
GBR_sample = samples[samples.population_code == 'GBR'].index
GBR_g = genotypes.subset(sel0=None, sel1=GBR_sample)
JPT_sample = samples[samples.population_code == 'JPT'].index
JPT_g = genotypes.subset(sel0=None, sel1=JPT_sample)
PJL_sample = samples[samples.population_code == 'PJL'].index
PJL_g = genotypes.subset(sel0=None, sel1=PJL_sample)
YRI_sample = samples[samples.population_code == 'YRI'].index
YRI_g = genotypes.subset(sel0=None, sel1=YRI_sample)
MXL_hom_ref = (MXL_g.count_hom_ref(axis=1))[:]
MXL_het = (MXL_g.count_het(axis=1))[:]
MXL_hom_alt = (MXL_g.count_hom_alt(axis=1))[:]
GBR_hom_ref = (GBR_g.count_hom_ref(axis=1))[:]
GBR_het = (GBR_g.count_het(axis=1))[:]
GBR_hom_alt = (GBR_g.count_hom_alt(axis=1))[:]
JPT_hom_ref = (JPT_g.count_hom_ref(axis=1))[:]
JPT_het = (JPT_g.count_het(axis=1))[:]
JPT_hom_alt = (JPT_g.count_hom_alt(axis=1))[:]
PJL_hom_ref = (PJL_g.count_hom_ref(axis=1))[:]
PJL_het = (PJL_g.count_het(axis=1))[:]
PJL_hom_alt = (PJL_g.count_hom_alt(axis=1))[:]
YRI_hom_ref = (YRI_g.count_hom_ref(axis=1))[:]
YRI_het = (YRI_g.count_het(axis=1))[:]
YRI_hom_alt = (YRI_g.count_hom_alt(axis=1))[:]


# Generating a dataframe for genotypes for each population
pos = variants['POS'][:]
MXL_genotypes = np.column_stack([pos,MXL_hom_ref,MXL_het,MXL_hom_alt])
MXL_genotypes = pd.DataFrame(MXL_genotypes)
GBR_genotypes = np.column_stack([pos,GBR_hom_ref,GBR_het,GBR_hom_alt])
GBR_genotypes = pd.DataFrame(GBR_genotypes)
JPT_genotypes = np.column_stack([pos,JPT_hom_ref,JPT_het,JPT_hom_alt])
JPT_genotypes = pd.DataFrame(JPT_genotypes)
PJL_genotypes = np.column_stack([pos,PJL_hom_ref,PJL_het,PJL_hom_alt])
PJL_genotypes = pd.DataFrame(PJL_genotypes)
YRI_genotypes = np.column_stack([pos,YRI_hom_ref,YRI_het,YRI_hom_alt])
YRI_genotypes = pd.DataFrame(YRI_genotypes)
columns = ['hom_ref','het','hom_alt']
MXL_genotypes.columns = ('POS','hom_ref','het','hom_alt')
MXL_genotypes['genotypes'] = MXL_genotypes[columns].to_dict(orient='records')
MXL_genotypes = MXL_genotypes.drop(columns=columns)
GBR_genotypes.columns = ('POS','hom_ref','het','hom_alt')
GBR_genotypes['genotypes'] = GBR_genotypes[columns].to_dict(orient='records')
GBR_genotypes = GBR_genotypes.drop(columns=columns)
JPT_genotypes.columns = ('POS','hom_ref','het','hom_alt')
JPT_genotypes['genotypes'] = JPT_genotypes[columns].to_dict(orient='records')
JPT_genotypes = JPT_genotypes.drop(columns=columns)
PJL_genotypes.columns = ('POS','hom_ref','het','hom_alt')
PJL_genotypes['genotypes'] = PJL_genotypes[columns].to_dict(orient='records')
PJL_genotypes = PJL_genotypes.drop(columns=columns)
YRI_genotypes.columns = ('POS','hom_ref','het','hom_alt')
YRI_genotypes['genotypes'] = YRI_genotypes[columns].to_dict(orient='records')
YRI_genotypes = YRI_genotypes.drop(columns=columns)


# Generating allele frequency dataframe for each population
columns = ['ref','alt']
MXL_ac = np.column_stack([id, pos, MXL_ac])
MXL_ac = pd.DataFrame(MXL_ac)
MXL_ac.columns = ('ID','POS','ref','alt')
MXL_ac['alleles'] = MXL_ac[columns].to_dict(orient='records')
MXL_ac = MXL_ac.drop(columns=columns)
GBR_ac = np.column_stack([id, pos,GBR_ac])
GBR_ac = pd.DataFrame(GBR_ac)
GBR_ac.columns = ('ID','POS','ref','alt')
GBR_ac['alleles'] = GBR_ac[columns].to_dict(orient='records')
GBR_ac = GBR_ac.drop(columns=columns)
JPT_ac = np.column_stack([id, pos, JPT_ac])
JPT_ac = pd.DataFrame(JPT_ac)
JPT_ac.columns = ('ID','POS','ref','alt')
JPT_ac['alleles'] = JPT_ac[columns].to_dict(orient='records')
JPT_ac = JPT_ac.drop(columns=columns)
PJL_ac = np.column_stack([id, pos, PJL_ac])
PJL_ac = pd.DataFrame(PJL_ac)
PJL_ac.columns = ('ID','POS','ref','alt')
PJL_ac['alleles'] = PJL_ac[columns].to_dict(orient='records')
PJL_ac = PJL_ac.drop(columns=columns)
YRI_ac = np.column_stack([id, pos, YRI_ac])
YRI_ac = pd.DataFrame(YRI_ac)
YRI_ac.columns = ('ID','POS','ref','alt')
YRI_ac['alleles'] = YRI_ac[columns].to_dict(orient='records')
YRI_ac = YRI_ac.drop(columns=columns)

# Final population dataframe to csv
MXL_data = MXL_ac
MXL_data['genotypes'] = MXL_genotypes['genotypes'].values
GBR_data = GBR_ac
GBR_data['genotypes'] = GBR_genotypes['genotypes'].values
JPT_data = JPT_ac
JPT_data['genotypes'] = JPT_genotypes['genotypes'].values
PJL_data = PJL_ac
PJL_data['genotypes'] = PJL_genotypes['genotypes'].values
YRI_data = YRI_ac
YRI_data['genotypes'] = YRI_genotypes['genotypes'].values



def genarray_to_list(array, data):
    lst = list(array)
    new_list = []
    for x in lst:
        y=x.tolist()
        new_list.append(y)
    string_list = []
    for x in new_list:
        x = str(x)
        string_list.append(x)
    data['genotype_array'] = string_list
    return data

YRI_data = genarray_to_list(YRI_g, YRI_data)
YRI_data.to_csv('YRI_data.tsv',sep='\t',index=False)
MXL_data = genarray_to_list(MXL_g, MXL_data)
MXL_data.to_csv('MXL_data.tsv',sep='\t',index=False)
JPT_data = genarray_to_list(JPT_g, JPT_data)
JPT_data.to_csv('JPT_data.tsv',sep='\t',index=False)
PJL_data = genarray_to_list(PJL_g, PJL_data)
PJL_data.to_csv('PJL_data.tsv',sep='\t',index=False)
GBR_data = genarray_to_list(GBR_g, GBR_data)
GBR_data.to_csv('GBR_data.tsv',sep='\t',index=False)

