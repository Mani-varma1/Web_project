from flask import render_template, request, Blueprint 
from VCF_website.query_results.routes import loading
from VCF_website.main.forms import SearchRs,SearchGene,SearchPos
# from VCF_website.results.

main = Blueprint("main",__name__)

@main.route("/")
@main.route("/home")
def home():
    return render_template('home.html')


@main.route("/about")
def about():
    return render_template('about.html', title='About')


@main.route("/search", methods=['GET', 'POST'])
def search():
    """ instantiating the forms seperately from the objects created in forms.py"""
    form1 = SearchPos()
    form2 = SearchRs()
    form3 = SearchGene()

    """ Check if form 1 was submmited and validations passsed"""
    if form1.submit.data and form1.validate_on_submit():
        """ Creating a Dictionary with form data if user searches via location
            for ease of identifying data on other pages.        
        """
        chromosome_position = {
            "chr": form1.select.data,
            "start_pos": form1.start_pos.data,
            "end_pos": form1.end_pos.data
        }        
        """ Instead of loading the route we can directly call in the function so the datatype is maintained
            else the dictionary turns into a string type.
        """
        return loading(search=chromosome_position)

    
        """
        Check if form 2 was submmited and validations passsed. Values are of string class,
        or a list of strings if multiple values seperated by commas are passed
        """
    elif form2.rs_search.data and form2.validate_on_submit():
        return loading(search=form2.rs_val.data)

        
        """
        Check if form 3 was submmited and validations passsed. Values are of string class,
        or a list of strings if multiple values seperated by commas are passed
        """
    elif form3.gene_search.data and form3.validate_on_submit():
        return loading(search=form3.gene.data)
    return render_template('search.html', title='About', form1=form1, form2=form2, form3=form3)




""" Redirects to documentation page""""
@main.route("/help", methods=['GET', 'POST'])
def help():
    return render_template('help.html', title='Contact')



""" Create Custom Error Pages for improving user experience if non existing URL
    was searched or if unexpected issue arise
"""
# Invalid URL
@main.errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404


# Internal Server Error
@main.errorhandler(500)
def page_not_found(e):
    return render_template("500.html"), 500
