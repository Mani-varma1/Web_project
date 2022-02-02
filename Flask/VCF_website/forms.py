from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField, IntegerField
from wtforms.validators import DataRequired, Length, Email, Optional
from VCF_website.models import VCF_db

class ContactForm(FlaskForm):
    name = StringField(label='Name', validators=[DataRequired()])
    email = StringField('Email', validators=[DataRequired(), Email()])
    message = StringField('Message',validators=[DataRequired(), Length(min=1, message=500)])
    submit = SubmitField('Send')


class SearchPos(FlaskForm):
    choices = [('CHR', 'CHR'),('1', '1'), ('22', '22')]
    select = SelectField('Select Chromosome:', choices=choices)
    start_pos = IntegerField('Starting position', validators=[DataRequired()])
    end_pos = IntegerField('Ending position', validators=[Optional()])
    submit = SubmitField('Search')
    
    
    
class SearchRs(FlaskForm):
    rs_val = StringField('rs value',validators=[DataRequired()])
    rs_search = SubmitField('Search')


class SearchGene(FlaskForm):
    gene = StringField('Gene Name',validators=[DataRequired()])
    gene_search = SubmitField('Search')
