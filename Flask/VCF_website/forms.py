from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField
from wtforms.validators import DataRequired, Length, Email
from VCF_website.models import VCF_db

class ContactForm(FlaskForm):
    name = StringField(label='Name', validators=[DataRequired()])
    email = StringField('Email', validators=[DataRequired(), Email()])
    message = StringField('Message',validators=[DataRequired(), Length(min=1, message=500)])
    submit = SubmitField('Sign Up')


class SearchForVCF(FlaskForm):
    choices = [('CHR', 'CHR'),('1', '1'), ('2', '2')]
    select = SelectField('Search for Kinase:', choices=choices)
    search = StringField('Name',validators=[DataRequired()])