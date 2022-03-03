import pandas as pd 
import numpy as np

altered = pd.read_csv('rs_id.txt')
altered
original_genes = pd.read_csv('rs_genes.txt',sep='\t',header=None)
original_genes
original_genes.columns = ['ID','gene_name']
original_none = pd.read_csv('rs_missing.txt')
original_none
original_none['gene_name'] = np.nan

total = pd.concat([original_genes,original_none],ignore_index=True)
total
altered.columns = ['ID']
final = pd.merge(altered,total, on='ID')
final.to_csv('final_gene_queries.tsv',sep='\t',index=False)
query = pd.read_table('query.tsv',sep='\t')
query = pd.DataFrame(query)
new_df = pd.merge(query,final, on='ID',how='left')
new_df
new_df.to_csv('query.tsv',sep='\t',index=False)