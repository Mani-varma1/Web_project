from flask_wtf import FlaskForm
from wtforms import SubmitField, SelectField, IntegerField, SelectMultipleField
from wtforms.validators import DataRequired, Optional, ValidationError
from markupsafe import Markup
from wtforms.widgets.core import html_params, ListWidget, CheckboxInput


class MultiCheckboxField(SelectMultipleField):
    """
    SelectMultipleField just creates a table of all availiable options to be selected.
    To improve user experience, we have combined it with the checkbox widgets from Wtforms
    allowing us to use check boxes, which is more intuative.
    """
    widget = ListWidget(prefix_label=False)
    option_widget =CheckboxInput()

    

class PopulationStatistics(FlaskForm):
    """
    Form provides users with checkboxs for selected the required populations and stats,
    with both fields requiring at least one selected option and a  custom validators 
    for checking if the user has selected more than 2 populations if FST stat was selected. 
    """
    choices_stats = [('Homozygosity','Homozygosity'),('Nucleotide Diversity','Nucleotide Diversity'),('Haplotype Diversity','Haplotype Diversity'),('Tajimas D','Tajimas D'),('FST','FST')]
    # render_kw provides css styles so there are no bullet points next to the checkbox as we are using unordered lists
    stats = MultiCheckboxField('Select Statistics:', choices=choices_stats, validators=[DataRequired(message="Please select atleast one statistic")],render_kw={'style': 'height: fit-content; list-style: none;'})
    choices_pop = [('GBR','British'),('JPT','Japanese'),('MXL','Mexican'),('PJL','Punjabi'),('YRI','Yoruba')]
    # render_kw provides css styles so there are no bullet points next to the checkbox as we are using unordered lists
    populations = MultiCheckboxField('Select Population:', choices=choices_pop, validators=[DataRequired(message="Please select a population")],render_kw={'style': 'height: fit-content; list-style: none;'})
    bin_size = IntegerField("Bin size:", validators=[DataRequired()])
    step_size = IntegerField("Step size (Optioanal):",validators=[Optional()])
    pop_stat = SubmitField("Calculate")
    
    # Custom validator to check if more than 2 populations are selected if fst was picked
    def validate_populations(self,populations):
        if ('FST' in  self.stats.data) and (len(populations.data)<2):
            raise ValidationError('Need atleast two populations for FST')
