from flask import render_template,url_for ,flash, redirect, request
from importlib_metadata import email
from VCF_website import app
from VCF_website.forms import ContactForm,SearchPos, SearchRs,SearchGene
from VCF_website.models import VCF_db





# @app.before_first_request
# def create_tables():
#     db.create_all()


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
        return redirect(url_for('results'))
    elif form2.rs_search.data and form2.validate_on_submit():
        return redirect(url_for('home'))
    elif form3.gene_search.data and form3.validate_on_submit():
        return redirect(url_for('results'))
    return render_template('search.html', title='About', form1=form1, form2=form2, form3=form3)



@app.route("/results",methods=['GET','POST'])
def results():
    return render_template('results.html', title='Results')
    

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


