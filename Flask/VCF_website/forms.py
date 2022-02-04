from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField, IntegerField
from wtforms.validators import DataRequired, Length, Email, Optional, ValidationError,InputRequired
from VCF_website.models import VCF_db
from markupsafe import Markup
from wtforms.widgets.core import html_params

class ContactForm(FlaskForm):
    name = StringField(label='Name', validators=[DataRequired()])
    email = StringField('Email', validators=[DataRequired(), Email()])
    message = StringField('Message',validators=[DataRequired(), Length(min=1, message=500)])
    submit = SubmitField('Send')


class CustomSelect:
    """
    Renders a select field allowing custom attributes for options.
    Expects the field to be an iterable object of Option fields.
    The render function accepts a dictionary of option ids ("{field_id}-{option_index}")
    which contain a dictionary of attributes to be passed to the option.

    Using this custom widget to disable CHR as an option so user cannot select
    """

    def __init__(self, multiple=False):
        self.multiple = multiple

    def __call__(self, field, option_attr=None, **kwargs):
        if option_attr is None:
            option_attr = {}
        kwargs.setdefault("id", field.id)
        if self.multiple:
            kwargs["multiple"] = True
        if "required" not in kwargs and "required" in getattr(field, "flags", []):
            kwargs["required"] = True
        html = ["<select %s>" % html_params(name=field.name, **kwargs)]
        for option in field:
            attr = option_attr.get(option.id, {})
            html.append(option(**attr))
        html.append("</select>")
        return Markup("".join(html))



class SearchPos(FlaskForm):
    choices = [("", "CHR"),('1', '1'), ('22', '22'),]
    select = SelectField("Select Chromosome", choices=choices, validators=[DataRequired(message="Please select a chromosome")],widget=CustomSelect(),default="",)
    start_pos = IntegerField("Starting position", validators=[DataRequired()])
    end_pos = IntegerField("Ending position", validators=[Optional()])
    submit = SubmitField("Search")

    
    
class SearchRs(FlaskForm):
    rs_val = StringField("rs value",validators=[DataRequired()])
    rs_search = SubmitField("Search")


class SearchGene(FlaskForm):
    gene = StringField("Gene Name",validators=[DataRequired()])
    gene_search = SubmitField("Search")
