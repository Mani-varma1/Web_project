import pandas as pd 
import sqlite3 as sql 
import numpy as np
import zipfile

with zipfile.ZipFile('data_for_populating.zip', 'r') as my_zip:
    my_zip.extractall('data')


conn = sql.connect('VCF_website/site.db')

c = conn.cursor()

mxl = pd.read_csv('data/MXL_data.tsv',sep='\t')
gbr = pd.read_csv('data/GBR_data.tsv',sep='\t')
jpt = pd.read_csv('data/JPT_data.tsv',sep='\t')
pjl = pd.read_csv('data/PJL_data.tsv',sep='\t')
yri = pd.read_csv('data/YRI_data.tsv',sep='\t')
query = pd.read_csv('data/query.tsv',sep='\t')

for row in query.itertuples():
    c.execute("""
            INSERT INTO query_search (chrom,rs_val,pos,gene_name,ref_allele,alt_allele)
            VALUES (?,?,?,?,?,?)
            """,(22,row.ID,row.POS,row.gene_name,row.REF,row.ALT)
    )

for row in mxl.itertuples():
    c.execute("""
            INSERT INTO snp_MXL (geno_freq,allele_freq,rs_val_id,genotypes)
            VALUES (?,?,?,?)
            """,(row.genotypes,row.alleles,row.ID,row.genotype_array)
    )

for row in gbr.itertuples():
    c.execute("""
            INSERT INTO snp_GBR (geno_freq,allele_freq,rs_val_id,genotypes)
            VALUES (?,?,?,?)
            """,(row.genotypes,row.alleles,row.ID,row.genotype_array)
    )

for row in jpt.itertuples():
    c.execute("""
            INSERT INTO snp_JPT (geno_freq,allele_freq,rs_val_id,genotypes)
            VALUES (?,?,?,?)
            """,(row.genotypes,row.alleles,row.ID,row.genotype_array)
    )

for row in pjl.itertuples():
    c.execute("""
            INSERT INTO snp_PJL (geno_freq,allele_freq,rs_val_id,genotypes)
            VALUES (?,?,?,?)
            """,(row.genotypes,row.alleles,row.ID,row.genotype_array)
    )
for row in yri.itertuples():
    c.execute("""
            INSERT INTO snp_YRI (geno_freq,allele_freq,rs_val_id,genotypes)
            VALUES (?,?,?,?)
            """,(row.genotypes,row.alleles,row.ID,row.genotype_array)
    )

conn.commit()

conn.close()
