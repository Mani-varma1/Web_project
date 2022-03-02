import pandas as pd 
import sqlite3 as sql 
import numpy as np

conn = sql.connect('site.db')

c = conn.cursor()

mxl = pd.read_csv('MXL_data.tsv',sep='\t')
gbr = pd.read_csv('GBR_data.tsv',sep='\t')
jpt = pd.read_csv('JPT_data.tsv',sep='\t')
pjl = pd.read_csv('PJL_data.tsv',sep='\t')
yri = pd.read_csv('YRI_data.tsv',sep='\t')
query = pd.read_csv('query.tsv',sep='\t')

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
