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
        """For Multi rsid or gene list, gene limited to just 5. if conditions not met
            redirects to search with a flashed message
        """
        if ',' in variable:
            try:
                
                variable = [i.strip() for i in variable.split(',')]
                rs_lst = []
                gene_lst=[]
                for i in variable:
                    if i.startswith('rs') and isinstance(int(i[2:]),int):
                        print(i)
                        print('hi')
                        rs_lst.append(i)
                    elif len(variable)<=5:
                        gene_lst.append(i)
                    elif len(variable)>5:
                        flash("Sorry Only 5 genes allowed", 'danger')
                        return redirect(url_for('search'))
                    else:
                        raise Exception
                if len(rs_lst)>0 and len(gene_lst)>0:
                    raise Exception
                else:
                    if len(rs_lst) != 0:
                        print(' DO query')
                        return '<h1>Temp RS Holder</h1>'
                    else:
                        print('Do query')
                        return '<h1>Temp GENE Holder</h1>'  
            except Exception:
                flash("Sorry please check your format and try again", 'danger')
                return redirect(url_for('search'))


        

        elif variable.startswith('rs') == True:
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

    for i in mxl:
        gf = ast.literal_eval(i['geno_freq'])
        gf_sum = sum(gf.values())
        gf_var = f"Hom-Ref:{round(gf['hom_ref']/gf_sum,2)}      Het:{round(gf['het']/gf_sum,2)}     Hom-Alt:{round(gf['hom_alt']/gf_sum,2)}"


        af = ast.literal_eval(i['allele_freq'])
        af_sum = sum(af.values())
        af_var = f"REF:{round(af['ref']/af_sum,2)}      ALT:{round(af['alt']/af_sum,2)}"

        i['geno_freq'] = gf_var
        i['allele_freq'] = af_var 

    form = PopulationStatistics()
    if request.method == "POST":
        if form.validate_on_submit():
            return redirect(url_for('stats',pops=form.populations.data, stats=form.stats.data))
    return render_template('results.html', Results=results, MXL=mxl, form=form)






## Comment this part









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







def win_pi_stats(positions,pop,bin_size,step_size=None):
    win_pi, pi_windows, nb, cts_pi = gstat.win_nuc_div(positions,pop,bin_size=bin_size,step_size=None)
    return win_pi, pi_windows, nb, cts_pi

def win_taj_stats(positions,pop,bin_size,step_size=None):
    win_tajima_D, win_taj, cts_taj = gstat.win_tajima_d(positions,pop,bin_size=bin_size,step_size=None)
    return win_tajima_D,win_taj,cts_taj

def moving_haplotype_div(pop,bin_size=100,step_size=None):
    moving_hap = gstat.moving_haplotype_div(pop=pop,bin_size=bin_size,step_size=step_size)
    return moving_hap



def win_fst_stats(positions,pop,bin_size,step_size=None):
    pass





@app.route("/stats/<pops>/<stats>")
def stats(pops, stats):
    results = json.loads(session['results'])
    gen_pos = [int(i['pos']) for i in results]

    first_col = results[0]
    last_col = results[-1]

    html_first_col = f"CHR:{first_col['chrom']} Start:{first_col['pos']} - End:{last_col['pos']}"
    html_gene = [set(i['gene_name']) for i in results if i['gene_name'] != None]


    """Summary Stats for each population"""
    """ GBR"""
    if 'GBR' in pops:
        gbr = json.loads(session['gbr'])
        gbr_gt_data, gbr_freq = decompress(gbr)

        """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""        
        gbr_homo, gbr_nuc_div, gbr_hap_div, gbr_taj_d = gstat.get_main_stats(pop=gbr_gt_data,freq_data =gbr_freq)
        gbr_stats = [gbr_homo, gbr_nuc_div, gbr_hap_div, gbr_taj_d]

        """Windowed Satats"""
        """PI"""
        gbr_win_pi, gbr_pi_windows, gbr_nb, gbr_cts_pi = win_pi_stats(positions=gen_pos,pop=gbr_gt_data,bin_size=10,step_size=None)

        """Tajimas D"""
        gbr_win_taj_D, gbr_win_taj, gbr_cts_taj= win_taj_stats(positions=gen_pos,pop=gbr_gt_data,bin_size=10,step_size=None)

        """ Haplotype"""
        gbr_mov_hap = moving_haplotype_div(pop=gbr_gt_data,bin_size=100,step_size=None)
    else:
        gbr = None

    """JPT"""
    if 'JPT' in pops:
        jpt = json.loads(session['jpt'])
        jpt_gt_data, jpt_freq = decompress(jpt)

        """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
        jpt_homo, jpt_nuc_div, jpt_hap_div, jpt_taj_d = gstat.get_main_stats(pop=jpt_gt_data,freq_data =jpt_freq)
        jpt_stats = [jpt_homo, jpt_nuc_div, jpt_hap_div, jpt_taj_d]

        """Windowed Satats"""

        """PI"""
        jpt_win_pi, jpt_pi_windows, jpt_nb, jpt_cts_pi= win_pi_stats(positions=gen_pos,pop=jpt_gt_data,bin_size=10,step_size=None)

        """Tajimas D"""
        jpt_win_taj_D, jpt_win_taj, jpt_cts_taj= win_taj_stats(positions=gen_pos,pop=jpt_gt_data,bin_size=10,step_size=None)

        """ Haplotype"""
        jpt_mov_hap = moving_haplotype_div(pop=jpt_gt_data,bin_size=100,step_size=None)
    else:
        jpt = None

    
    
    """MXL"""
    if 'MXL' in pops:
        mxl = json.loads(session['mxl'])
        mxl_gt_data, mxl_freq = decompress(mxl)
        """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
        mxl_homo, mxl_nuc_div, mxl_hap_div, mxl_taj_d = gstat.get_main_stats(pop=mxl_gt_data,freq_data =mxl_freq)
        mxl_stats = [mxl_homo, mxl_nuc_div, mxl_hap_div, mxl_taj_d]

        """Windowed Satats"""
        """PI"""
        mxl_win_pi, mxl_pi_windows, mxl_nb, mxl_cts_pi= win_pi_stats(positions=gen_pos,pop=mxl_gt_data,bin_size=10,step_size=None)

        """Tajimas D"""
        mxl_win_taj_D, mxl_win_taj, mxl_cts_taj= win_taj_stats(positions=gen_pos,pop=mxl_gt_data,bin_size=10,step_size=None)

        """ Haplotype"""
        mxl_mov_hap = moving_haplotype_div(pop=mxl_gt_data,bin_size=100,step_size=None)
    else:
        mxl = None




    if 'PJL' in pops:
        pjl = json.loads(session['pjl'])
        pjl_gt_data, pjl_freq = decompress(pjl)
        """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
        pjl_homo, pjl_nuc_div, pjl_hap_div, pjl_taj_d = gstat.get_main_stats(pop=pjl_gt_data,freq_data =pjl_freq)
        pjl_stats = [pjl_homo, pjl_nuc_div, pjl_hap_div, pjl_taj_d]

        """Windowed Satats"""
        """PI"""
        pjl_win_pi, pjl_pi_windows, pjl_nb, pjl_cts_pi= win_pi_stats(positions=gen_pos,pop=pjl_gt_data,bin_size=10,step_size=None)

        """Tajimas D"""
        pjl_win_taj_D, pjl_win_taj, pjl_cts_taj= win_taj_stats(positions=gen_pos,pop=pjl_gt_data,bin_size=10,step_size=None)

        """ Haplotype"""
        pjl_mov_hap = moving_haplotype_div(pop=pjl_gt_data,bin_size=100,step_size=None)
    else:
        pjl = None



    if 'YRI' in pops:
        yri = json.loads(session['yri'])
        yri_gt_data, yri_freq = decompress(yri)
        """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
        yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d = gstat.get_main_stats(pop=yri_gt_data,freq_data =yri_freq)
        yri_stats = [yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d]

        """Windowed Satats"""
        """PI"""
        yri_win_pi, yri_pi_windows, yri_nb, yri_cts_pi= win_pi_stats(positions=gen_pos,pop=yri_gt_data,bin_size=10,step_size=None)

        """Tajimas D"""
        yri_win_taj_D, yri_win_taj, yri_cts_taj= win_taj_stats(positions=gen_pos,pop=yri_gt_data,bin_size=10,step_size=None)

        """ Haplotype"""
        yri_mov_hap = moving_haplotype_div(pop=yri_gt_data,bin_size=100,step_size=None)
    else:
        yri = None



    gt_dict = {}
    gt_freq = {}

    if gbr:
        gt_dict['GBR'] = gbr_gt_data
        gt_freq['GBR'] = gbr_freq
    
    if jpt:
        gt_dict['JPT'] = jpt_gt_data
        gt_freq['JPT'] = jpt_freq

    if mxl:
        gt_dict['MXL'] = mxl_gt_data
        gt_freq['MXL'] = mxl_freq
    
    if pjl:
        gt_dict['PJL'] =pjl_gt_data
        gt_freq['PJL'] = pjl_freq

    if yri:
        gt_dict['YRI'] = pjl_gt_data
        gt_freq['YRI'] = yri_freq

    


    all_fstat = gstat.get_fstat(paris=pops, gt_dict = gt_dict)
    all_pops =[gt_dict[i] for i in pops]
    
    print(all_pops)




    return render_template('stats.html', stats=stats, populations=pops, results=results, )







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

