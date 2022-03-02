from VCF_website import db
# from VCF_website import ma


class query_search(db.Model):
    chrom = db.Column(db.String(20), unique=False, nullable=False)
    rs_val = db.Column(db.String(20), primary_key=True, unique=True, nullable=False)
    pos = db.Column(db.String(120), unique=True, nullable=False)
    gene_name = db.Column(db.String(120), unique=False)
    ref_allele = db.Column(db.String(20), unique=False, nullable=False)
    alt_allele = db.Column(db.String(20), unique=False, nullable=False)
    mxl = db.relationship('snp_MXL', backref='rs_val', lazy='select')
    gbr = db.relationship('snp_GBR', backref='rs_val', lazy='select')
    jpt = db.relationship('snp_JPT', backref='rs_val', lazy='select')
    pjl = db.relationship('snp_PJL', backref='rs_val', lazy='select')
    yri = db.relationship('snp_YRI', backref='rs_val', lazy='select')

    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns if c.name != 'id'}

class snp_MXL(db.Model):
    
    id = db.Column(db.Integer, primary_key=True)
    geno_freq = db.Column(db.String(120), unique=False, nullable=False)
    allele_freq = db.Column(db.String(120), unique=False, nullable=False)
    rs_val_id = db.Column(db.String, db.ForeignKey('query_search.rs_val'))
    genotypes = db.Column(db.String, unique=False)
    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns if c.name != 'id'}


class snp_GBR(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    geno_freq = db.Column(db.String(120), unique=False, nullable=False)
    allele_freq = db.Column(db.String(120), unique=False, nullable=False)
    rs_val_id = db.Column(db.String, db.ForeignKey('query_search.rs_val'))
    genotypes = db.Column(db.String, unique=False)
    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns if c.name != 'id'}


class snp_PJL(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    geno_freq = db.Column(db.String(120), unique=False, nullable=False)
    allele_freq = db.Column(db.String(120), unique=False, nullable=False)
    rs_val_id = db.Column(db.String, db.ForeignKey('query_search.rs_val'))
    genotypes = db.Column(db.String, unique=False)

    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns if c.name != 'id'}


class snp_JPT(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    geno_freq = db.Column(db.String(120), unique=False, nullable=False)
    allele_freq = db.Column(db.String(120), unique=False, nullable=False)
    rs_val_id = db.Column(db.String, db.ForeignKey('query_search.rs_val'))
    genotypes = db.Column(db.String, unique=False)
    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns if c.name != 'id'}


class snp_YRI(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    geno_freq = db.Column(db.String(120), unique=False, nullable=False)
    allele_freq = db.Column(db.String(120), unique=False, nullable=False)
    rs_val_id = db.Column(db.String, db.ForeignKey('query_search.rs_val'))
    genotypes = db.Column(db.String, unique=False)
    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns if c.name != 'id'}
