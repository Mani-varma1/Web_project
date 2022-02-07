from VCF_website import db



class query_search(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column(db.String(20), unique=False, nullable=False)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    pos = db.Column(db.String(120), unique=True, nullable=False)
    gene_name = db.Column(db.String(120), unique=False)
    ref_allele = db.Column(db.String(20), unique=False, nullable=False)
    alt_allele = db.Column(db.String(20), unique=False, nullable=False)

    def __repr__(self):
        return 'chrom: {0} | RSID: {1} | Position: {2} | Reference: {3} | Alternate: {4} |'.format(self.chrom,self.rs_val,self.pos,self.ref_allele,self.alt_allele)


class sample_info(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    sample_id = db.Column(db.String(20), nullable=False)
    population_code = db.Column(db.String(20), nullable=False, unique=False)
    superpopulation_code = db.Column(db.String(20), nullable=False, unique=False)

class genotype_sample(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    sample_id = db.Column(db.String(20), nullable=False, unique=False)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    genotype = db.Column(db.String(20), nullable=False, unique=False)

class snp_MXL(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    pos = db.Column(db.String(120), unique=True, nullable=False)
    ref_allele = db.Column(db.String(20), unique=False, nullable=False)
    alt_allele = db.Column(db.String(20), unique=False, nullable=False)
    geno_freq = db.Column(db.String(20), unique=False, nullable=False)
    allele_freq = db.Column(db.String(20), unique=False, nullable=False)

class snp_GBR(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    pos = db.Column(db.String(120), unique=True, nullable=False)
    ref_allele = db.Column(db.String(20), unique=False, nullable=False)
    alt_allele = db.Column(db.String(20), unique=False, nullable=False)
    geno_freq = db.Column(db.String(20), unique=False, nullable=False)
    allele_freq = db.Column(db.String(20), unique=False, nullable=False)

class snp_PJL(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    pos = db.Column(db.String(120), unique=True, nullable=False)
    ref_allele = db.Column(db.String(20), unique=False, nullable=False)
    alt_allele = db.Column(db.String(20), unique=False, nullable=False)
    geno_freq = db.Column(db.String(20), unique=False, nullable=False)
    allele_freq = db.Column(db.String(20), unique=False, nullable=False)

class snp_JPT(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    pos = db.Column(db.String(120), unique=True, nullable=False)
    ref_allele = db.Column(db.String(20), unique=False, nullable=False)
    alt_allele = db.Column(db.String(20), unique=False, nullable=False)
    geno_freq = db.Column(db.String(20), unique=False, nullable=False)
    allele_freq = db.Column(db.String(20), unique=False, nullable=False)

class snp_YRI(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    pos = db.Column(db.String(120), unique=True, nullable=False)
    ref_allele = db.Column(db.String(20), unique=False, nullable=False)
    alt_allele = db.Column(db.String(20), unique=False, nullable=False)
    geno_freq = db.Column(db.String(20), unique=False, nullable=False)
    allele_freq = db.Column(db.String(20), unique=False, nullable=False)
    
