from VCF_website.models import snp_GBR,snp_JPT,snp_MXL,snp_PJL,snp_YRI,query_search
from flask import session
import json
import ast


''' This function is used to run queries for each population, depending on what search parameters were used. 
Uses SQLAlchemy to structure the queries.
Query results are then stored in sessions as JSON objects which can then be loaded to the results and stats pages.
'''
def pop_data(results,variable):
    mxl = []
    gbr = []
    jpt = []
    pjl = []
    yri = []
    if isinstance(variable,dict):
        if variable['end_pos'] != None:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
        else:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
    elif isinstance(variable,list):
            if variable[0].startswith('rs') and isinstance(int((variable[0])[2:]),int):
                mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all()
                gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
                jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
                pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
                yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
            else:
                mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all()
                gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
                jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
                pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
                yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
    elif isinstance(variable,str):
        if variable.startswith('rs') == True:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all() 
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all() 
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all() 
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all()
        else:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all() 
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all() 
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all() 
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all()

    session['results'] = json.dumps([i.to_dict() for i in results])
    session['mxl'] = json.dumps([i.to_dict() for i in mxl])
    session['gbr'] = json.dumps([i.to_dict() for i in gbr])
    session['jpt'] = json.dumps([i.to_dict() for i in jpt])
    session['pjl'] = json.dumps([i.to_dict() for i in pjl])
    session['yri'] = json.dumps([i.to_dict() for i in yri])
    
    return None


''' This function converts the genotype and allele counts to frequencies to be displayed on the results page. '''
def convert_freq(pop):
    for i in pop:
        gf = ast.literal_eval(i['geno_freq'])
        gf_sum = sum(gf.values())
        gf_var = f"Hom-Ref:{round(gf['hom_ref']/gf_sum,2)}         Het:{round(gf['het']/gf_sum,2)}        Hom-Alt:{round(gf['hom_alt']/gf_sum,2)}"


        af = ast.literal_eval(i['allele_freq'])
        af_sum = sum(af.values())
        af_var = f"REF:{round(af['ref']/af_sum,2)}      ALT:{round(af['alt']/af_sum,2)}"

        i['geno_freq'] = gf_var
        i['allele_freq'] = af_var 

    return None
