from flask_wtf import FlaskForm
from wtforms import SubmitField, SelectField, IntegerField, SelectMultipleField
from wtforms.validators import DataRequired, Optional, ValidationError
from markupsafe import Markup
from wtforms.widgets.core import html_params, ListWidget, CheckboxInput

"""
SelectMultipleField just creates a table of all availiable options to be selected.
To improve user experience, we have combined it with the checkbox widgets from Wtforms
allowing us to use check boxes, which is more intuative.
"""
class MultiCheckboxField(SelectMultipleField):
    widget = ListWidget(prefix_label=False)
    option_widget =CheckboxInput()

    
"""
Form provides users with checkboxs for selected the required populations and stats,
with both fields requiring at least one selected option and a  custom validators 
for checking if the user has selected more than 2 populations if FST stat was selected. 
"""
class PopulationStatistics(FlaskForm):
    choices_stats = [('Homozygosity','Homozygosity'),('Nucleotide Diversity','Nucleotide Diversity'),('Haplotype Diversity','Haplotype Diversity'),('Tajimas D','Tajimas D'),('FST','FST')]
    stats = MultiCheckboxField('Select Statistics:', choices=choices_stats, validators=[DataRequired(message="Please select atleast one statistic")])
    choices_pop = [('GBR','GBR'),('JPT','JPT'),('MXL','MXL'),('PJL','PJL'),('YRI','YRI')]
    populations = MultiCheckboxField('Select Population:', choices=choices_pop, validators=[DataRequired(message="Please select a population")])
    bin_size = IntegerField("Bin size:", validators=[DataRequired()])
    step_size = IntegerField("Step size (Optioanal):",validators=[Optional()])
    pop_stat = SubmitField("Calculate")

    def validate_populations(self,populations):
        if ('FST' in  self.stats.data) and (len(populations.data)<2):
            raise ValidationError('Need atleast two populations for FST')
