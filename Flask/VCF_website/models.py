from VCF_website import db


class VCF_db(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    Chrom = db.Column(db.String(20), unique=True, nullable=False)
    rs_val = db.Column(db.String(20), unique=True, nullable=False)
    Pos = db.Column(db.String(120), unique=True, nullable=False)
    qual = db.Column(db.String(20), nullable=False, default='default.jpg')

    def __repr__(self):
        return f"User('{self.id}', '{self.Chrom}', '{self.Pos}','{self.rs_val}')"

