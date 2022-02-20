import code
import time
from flask import render_template, url_for, flash, redirect, request, session, make_response
from io import StringIO
from numpy import append
from werkzeug.wrappers import Response
from VCF_website import app,sess
from VCF_website.forms import ContactForm, SearchPos, SearchRs, SearchGene, PopulationStatistics
from VCF_website.models import query_search, snp_MXL, snp_GBR, snp_JPT, snp_PJL, snp_YRI
import VCF_website.genome_stats as gstat
import ast
import csv
import json

# @app.before_first_request
# def create_tables():
# db.create_all()


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






def pop_data(results,variable):
    mxl = []
    gbr = []
    jpt = []
    pjl = []
    yri = []
    if isinstance(variable,dict):
        if variable['end_pos'] != None:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            print('hi')
        else:
            for x in results:
                mxl = x.mxl
                gbr = x.gbr
                jpt = x.jpt
                pjl = x.pjl
                yri = x.yri
    else:
        for x in results:
            mxl = x.mxl
            gbr = x.gbr
            jpt = x.jpt
            pjl = x.pjl
            yri = x.yri


    session['results'] = json.dumps([i.to_dict() for i in results])
    session['mxl'] = json.dumps([i.to_dict() for i in mxl])
    session['gbr'] = json.dumps([i.to_dict() for i in gbr])
    session['jpt'] = json.dumps([i.to_dict() for i in jpt])
    session['pjl'] = json.dumps([i.to_dict() for i in pjl])
    session['yri'] = json.dumps([i.to_dict() for i in yri])
    
    return None








@app.route("/loading", methods=['GET', 'POST'])
def loading(search):
    variable = search
    session.clear()
    if isinstance(variable, dict):
        if variable["end_pos"] == None:
            results = query_search.query.filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            
            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))
            
            pop_data(results,variable)
            
            return redirect(url_for('results', title='Results', Results=results))

        else:
            results = query_search.query.filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()

            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))


            start_time = time.time()
 
            pop_data(results,variable)

            print("--- %s seconds pop_data_func() ---" % (time.time() - start_time))

            return redirect(url_for('results', title='Results'))
    else:
        if variable.startswith('rs') == True:
            results = query_search.query.filter(query_search.rs_val.like(variable)).all() 

            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))

            pop_data(results,variable)
            return redirect(url_for('results', title='Results'))
        
        else:
            results = query_search.query.filter(query_search.gene_name.like(variable)).all()

            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))

            
            pop_data(results,variable)
            
            return redirect(url_for('results', title='Results', Results=results))







@app.route("/results", methods=['GET', 'POST'])
def results():
    results = json.loads(session['results'])
    mxl = json.loads(session['mxl'])
    form = PopulationStatistics()
    if request.method == "POST":
        if form.validate_on_submit():
            return stats(pop_sel=form.populations.data, stats_sel=form.stats.data)
    return render_template('results.html', Results=results, MXL=mxl, form=form)






## Comment this part







def get_main_stats(pop,freq_data):
    homo = gstat.Homozygosity(freq_data)
    nuc_div = gstat.nuc_div(pop)
    taj_d = gstat.tajima_d(pop)
    hap_div = gstat.haplotype_div(pop)
    return homo,nuc_div,hap_div,taj_d

def get_fstat():
    pass


def win_stats():
    pass




def decompress(gt_arr):
    freq_data = []
    gt_data = []
    decomp_dict = {'a':'[0, 0]','b':'[0, 1]','c':'[1, 0]','d':'[1,1]'}
    for item in gt_arr:
        gt_arr_data = ast.literal_eval(item['genotypes'])
        freq_data.append(ast.literal_eval(item['geno_freq']))

        snp_data = ""
        for sample in gt_arr_data:
            rep = int(sample[0])
            val = decomp_dict[sample[1]]+','
            val = rep*val
            snp_data += val

        gt_data.append(list(ast.literal_eval(snp_data[:-1])))

    return gt_data, freq_data




@app.route("/stats")
def stats(pop_sel, stats_sel):
    pop_sel = pop_sel
    stats_sel = stats_sel
    results = json.loads(session['results'])
    gen_pos = [int(i['pos']) for i in results]



    first_col = results[0]
    last_col = results[-1]

    html_first_col = f"CHR:{first_col['chrom']} Start:{first_col['pos']} - End:{last_col['pos']}"
    html_gene = [set(i['gene_name']) for i in results if i['gene_name'] != None]

    print(html_first_col)
    print(html_gene)







    if 'GBR' in pop_sel:
        gbr = json.loads(session['gbr'])
        gbr_gt_data, gbr_freq = decompress(gbr)
        
        gbr_homo,gbr_nuc_div,gbr_hap_div,gbr_taj_d = get_main_stats(pop=gbr_gt_data,freq_data =gbr_freq)
        gbr_stats = [gbr_homo,gbr_nuc_div,gbr_hap_div,gbr_taj_d]

    else:
        gbr = None


    if 'JPT' in pop_sel:
        jpt = json.loads(session['jpt'])
        jpt_gt_data, jpt_freq = decompress(jpt)
        jpt_homo,jpt_nuc_div,jpt_hap_div,jpt_taj_d = get_main_stats(pop=jpt_gt_data,freq_data =jpt_freq)
        jpt_stats = [jpt_homo,jpt_nuc_div,jpt_hap_div,jpt_taj_d]
    else:
        jpt = None

    
    
    
    if 'MXL' in pop_sel:
        mxl = json.loads(session['mxl'])
        mxl_gt_data, mxl_freq = decompress(mxl)
        mxl_homo,mxl_nuc_div,mxl_hap_div,mxl_taj_d = get_main_stats(pop=mxl_gt_data,freq_data =mxl_freq)
        mxl_stats = [mxl_homo,mxl_nuc_div,mxl_hap_div,mxl_taj_d]
    else:
        mxl = None




    if 'PJL' in pop_sel:
        pjl = json.loads(session['pjl'])
        pjl_gt_data, pjl_freq = decompress(pjl)
        pjl_homo,pjl_nuc_div,pjl_hap_div,pjl_taj_d = get_main_stats(pop=pjl_gt_data,freq_data =pjl_freq)
        pjl_stats = [pjl_homo,pjl_nuc_div,pjl_hap_div,pjl_taj_d]
    else:
        pjl = None



    if 'YRI' in pop_sel:
        yri = json.loads(session['yri'])
        yri_gt_data, yri_freq = decompress(yri)
        yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d = get_main_stats(pop=yri_gt_data,freq_data =yri_freq)
        yri_stats = [yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d]
    else:
        yri = None



    overall_stat = []

    if gbr:
        overall_stat.append(gbr_gt_data())
    
    if jpt:
        overall_stat.append(jpt_gt_data())

    if mxl:
        temp_mxl = mxl_gt_data
    
    if pjl:
        overall_stat.append(pjl_gt_data)

    if yri:
        overall_stat.append(yri_gt_data)

    # for i in 





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
