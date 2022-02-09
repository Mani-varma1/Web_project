from re import S
from flask import render_template,url_for ,flash, redirect, request, session
from importlib_metadata import email
from VCF_website import app
from VCF_website.forms import ContactForm,SearchPos, SearchRs,SearchGene
from VCF_website.models import query_search
import ast




#@app.before_first_request
#def create_tables():
    #db.create_all()


@app.route("/")
@app.route("/home")
def home():
    return render_template('home.html')

@app.route("/about")
def about():
    return render_template('about.html', title='About Us')


@app.route("/search",methods=['GET','POST'])
def search():
    form1 = SearchPos()
    form2 = SearchRs()
    form3 = SearchGene()

    if form1.submit.data and form1.validate_on_submit():
        chromosome_position =  {
            "chr" : form1.select.data,
            "start_pos" : form1.start_pos.data,
            "end_pos" : form1.end_pos.data
        }
        return redirect(url_for('results',search=chromosome_position))
    elif form2.rs_search.data and form2.validate_on_submit():
        form2_rs = f"'{form2.rs_val.data}'"
        return redirect(url_for('results',search=form2_rs))
    elif form3.gene_search.data and form3.validate_on_submit():
        form3_gene = f"'{form3.gene.data}'"
        return redirect(url_for('results',search=form3_gene))
    return render_template('search.html', title='About', form1=form1, form2=form2, form3=form3)

def get_results(variable):
    variable = ast.literal_eval(variable)
    if isinstance(variable, dict):
        if variable["end_pos"] == None:
            results = query_search.query.filter(query_search.pos.like(int(variable["start_pos"]))).all()
            if len(results) == 0:
                results = "None"
            return render_template('results.html', title='Results', Results=results)
        else:
            results = query_search.query.filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).all()
            if len(results) == 0:
                results = "None"
            return render_template('results.html', title='Results', Results=results)
    else:
        if variable.startswith('rs') == True:
            results = query_search.query.filter(query_search.rs_val.like(variable)).all()
            if len(results) == 0:
                results = "None"
            return render_template('results.html', title='Results', Results=results)
        else: 
            results = query_search.query.filter(query_search.gene_name.like(variable)).all()
            if len(results) == 0:
                results = "None"
            return render_template('results.html', title='Results', Results=results)        
                

@app.route("/loading/<variable>",methods=["GET","POST"])
def loading(variable):
    if request.method == "GET":
        results = get_result(variable)
        print(type(results))
        #if len(results) == 0:
        #    results = "None"
        return redirect(url_for("results",variable=results))
    return render_template('loading.html', title='Loading')

@app.route("/results/<search>",methods=['GET','POST'])
def results(search):
    if request.method == "GET":
        variable = search
        variable = ast.literal_eval(variable)
        if isinstance(variable, dict):
            if variable["end_pos"] == None:
                results = query_search.query.filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
                return render_template('results.html', title='Results', Results=results)
            else:
                results = query_search.query.filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
                return render_template('results.html', title='Results', Results=results)
        else:
            if variable.startswith('rs') == True:
                results = query_search.query.filter(query_search.rs_val.like(variable)).all()
                return render_template('results.html', title='Results', Results=results)
            else: 
                results = query_search.query.filter(query_search.gene_name.like(variable)).all()
                return render_template('results.html', title='Results', Results=results) 

    

@app.route("/contact", methods=['GET','POST'])
def contact():
    form = ContactForm()
    if form.validate_on_submit():
            flash(f'You Query has been submmited', 'success')
            return redirect(url_for('home'))
    return render_template('contact.html', title='Contact', form=form)




# def save_picture(form_picture):
#     random_hex = secrets.token_hex(8)
#     _, f_ext = os.path.splitext(form_picture.filename)
#     picture_fn = random_hex+f_ext
#     picture_path = os.path.join(app.root_path, 'static/profile_pics', picture_fn)
#     form_picture.save(picture_path)
#     return picture_fn


