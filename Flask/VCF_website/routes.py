from flask import render_template, url_for, flash, redirect, request, session, make_response
from io import StringIO
from werkzeug.wrappers import Response
from VCF_website import app,sess,db
from VCF_website.forms import ContactForm, SearchPos, SearchRs, SearchGene, PopulationStatistics
from VCF_website.models import query_search, snp_MXL, snp_GBR, snp_JPT, snp_PJL, snp_YRI
import ast
import csv
import allel
import numpy as np
import json
import VCF_website.genome_stats as gstat
import time

@app.before_first_request
def create_tables():
    db.create_all()


@app.route("/")
@app.route("/home")
def home():
    return render_template('home.html')





@app.route("/about")
def about():
    return render_template('about.html', title='About Us')






@app.route("/search", methods=['GET', 'POST'])
def search():
    form1 = SearchPos()
    form2 = SearchRs()
    form3 = SearchGene()
    if form1.submit.data and form1.validate_on_submit():
        chromosome_position = {
            "chr": form1.select.data,
            "start_pos": form1.start_pos.data,
            "end_pos": form1.end_pos.data
        }
        return loading(search=chromosome_position)
    elif form2.rs_search.data and form2.validate_on_submit():
        return loading(search=form2.rs_val.data)
    elif form3.gene_search.data and form3.validate_on_submit():
        return loading(search=form3.gene.data)
    return render_template('search.html', title='About', form1=form1, form2=form2, form3=form3)






def pop_data(variable, results, *end_pos):
    mxl = []
    gbr = []
    jpt = []
    pjl = []
    yri = []
    if end_pos != None:
        # for x in results:
        mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
        gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
        jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
        pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
        yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
    else:
        for x in results:
            mxl = x.mxl
            gbr = x.gbr
            jpt = x.jpt
            pjl = x.pjl
            yri = x.yri
    return mxl, gbr, jpt, pjl, yri








@app.route("/loading", methods=['GET', 'POST'])
def loading(search):
    variable = search
    session.clear()
    if isinstance(variable, dict):
        if variable["end_pos"] == None:
            results = snp_MXL.query.filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            mxl, gbr, jpt, pjl, yri = pop_data(results, variable["end_pos"])
            return redirect(url_for('results', title='Results', Results=results))
        else:
            start_time = time.time() 
            results = query_search.query.filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            print("--- %s seconds ---" % (time.time() - start_time))
            print(results)
            mxl, gbr, jpt, pjl, yri = pop_data(variable, results, variable["end_pos"])
            print(gbr)
            print("--- %s seconds ---" % (time.time() - start_time))            
            session['results'] = json.dumps([i.to_dict() for i in results])
            session['mxl'] = json.dumps([i.to_dict() for i in mxl])
            session['gbr'] = json.dumps([i.to_dict() for i in gbr])
            session['jpt'] = json.dumps([i.to_dict() for i in jpt])
            session['pjl'] = json.dumps([i.to_dict() for i in pjl])
            session['yri'] = json.dumps([i.to_dict() for i in yri])
            print("--- %s seconds ---" % (time.time() - start_time))           
            return redirect(url_for('results', title='Results'))
    else:
        if variable.startswith('rs') == True:
            results = query_search.query.filter(query_search.rs_val.like(variable)).all()
            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))
            session['results'] = json.dumps([i.to_dict() for i in results])
            session['mxl'] = json.dumps([i.to_dict() for i in mxl])
            session['gbr'] = json.dumps([i.to_dict() for i in gbr])
            session['jpt'] = json.dumps([i.to_dict() for i in jpt])
            session['pjl'] = json.dumps([i.to_dict() for i in pjl])
            session['yri'] = json.dumps([i.to_dict() for i in yri])
            return redirect(url_for('results', title='Results'))
        else:
            results = query_search.query.filter(query_search.gene_name.like(variable)).all()
            mxl, gbr, jpt, pjl, yri = pop_data(results)
            return redirect(url_for('results', title='Results', Results=results, MXL=mxl, GBR=gbr, JPT=jpt, PJL=pjl, YRI=yri))







@app.route("/results", methods=['GET', 'POST'])
def results():
    results = json.loads(session['results'])
    mxl = json.loads(session['mxl'])
    # pop_array = []
    # temp = []
    # for x in mxl:
    #     temp.append(x['genotypes'])

    # print(temp)
    # gen_arr = allel.GenotypeArray(np.array(pop_array))
    form = PopulationStatistics()
    if request.method == "POST":
        if form.validate_on_submit():
            return stats(pop_sel=form.populations.data, stats_sel=form.stats.data)
    return render_template('results.html', Results=results, MXL=mxl, form=form)


## Comment this part


@app.route("/stats")
def stats(pop_sel, stats_sel):
    results = json.loads(session['results'])
    mxl = json.loads(session['mxl'])
    mxl_positions = [int(i['pos']) for i in results]
    # pop_array = []
    # # for x in mxl:
    # #     pop_array.append(ast.literal_eval(x['genotypes']))
    # # gen_arr = allel.GenotypeArray(np.array(pop_array))
    # # print(gen_arr)
    mxl_homo = gstat.Homozygosity(mxl)


    mxl_haplotype_div = gstat.haplotype_div(mxl,10)


    mxl_td_td,mxl_td_win,mxl_td_cts = gstat.tajima_d(positions= mxl_positions,pop= mxl)
    
    
    mxl_nu_di_pi,mxl_nu_di_win,mxl_nu_di_nb, mxl_nu_di_cts = gstat.nucleotide_div(positions=mxl_positions, pop=mxl)


    print(mxl_homo)
    print(mxl_haplotype_div)
    print(mxl_td_td,mxl_td_win,mxl_td_cts)
    print(mxl_nu_di_pi,mxl_nu_di_win,mxl_nu_di_nb, mxl_nu_di_cts)


    return render_template('stats.html', stats=stats_sel, populations=pop_sel, results=results, )


"""
Temp Download page with dummy data

"""
test_down = [{'chrom': '22', 'rs_val': 'rs587698813', 'pos': '16051164',
              'gene_name': None, 'ref_allele': 'G', 'alt_allele': 'A'}]


@app.route('/download')
def download():
    si = StringIO()
    fields = [
        'chrom',
        'rs_val',
        'pos',
        'gene_name',
        'ref_allele',
        'alt_allele'
    ]
    cw = csv.DictWriter(si, fieldnames=fields)
    cw.writeheader()
    for stats in test_down:
        cw.writerow(stats)
    output = make_response(si.getvalue())
    output.headers["Content-Disposition"] = "attachment; filename=stats.csv"
    output.headers["Content-type"] = "text/csv"
    return output


@app.route("/contact", methods=['GET', 'POST'])
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
