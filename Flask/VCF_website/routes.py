
from flask import render_template, url_for, flash, redirect, request, session, make_response,send_file
from io import StringIO,BytesIO
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
    return render_template('about.html', title='About')




@app.route("/search", methods=['GET', 'POST'])
def search():
    """ instantiating the forms seperately from the objects created in forms.py"""
    form1 = SearchPos()
    form2 = SearchRs()
    form3 = SearchGene()

    """ Check if form 1 was submmited and validations passsed"""
    if form1.submit.data and form1.validate_on_submit():
        chromosome_position = {
            "chr": form1.select.data,
            "start_pos": form1.start_pos.data,
            "end_pos": form1.end_pos.data
        }
        return loading(search=chromosome_position)

        """
        Check if form 2 was submmited and validations passsed
        """
    elif form2.rs_search.data and form2.validate_on_submit():
        return loading(search=form2.rs_val.data)

        
        """
        Check if form 3 was submmited and validations passsed
        """
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
        else:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()
    elif isinstance(variable,list):
            if variable[0].startswith('rs') and isinstance(int((variable[0])[2:]),int):
                mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all()
                gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
                jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
                pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
                yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.rs_val.in_(variable)).all() 
            else:
                mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all()
                gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
                jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
                pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
                yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.gene_name.in_(variable)).all() 
    elif isinstance(variable,str):
        if variable.startswith('rs') == True:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all() 
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all() 
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all() 
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.rs_val.like(variable)).all()
        else:
            mxl = snp_MXL.query.filter(snp_MXL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all()
            gbr = snp_GBR.query.filter(snp_GBR.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all() 
            jpt = snp_JPT.query.filter(snp_JPT.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all() 
            pjl = snp_PJL.query.filter(snp_PJL.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all() 
            yri = snp_YRI.query.filter(snp_YRI.rs_val_id == query_search.rs_val).filter(query_search.gene_name.like(variable)).all()

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
    """ Explicitly clear sessions"""
    session.clear()
    
    
    """ Check if the user input was using Chromosome locations"""

    """ For single location"""
    if isinstance(variable, dict):
        if variable["end_pos"] == None:
            results = query_search.query.filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()

            """ If nothing found redirect"""
            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))

            """ Run the function that parses the data from query object to python list object and assign them to sessions"""
            pop_data(results,variable)

            """ Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation"""
            flash("Please read the documentation for appropriate parameters", 'info')            
            return redirect(url_for('results', title='Results', Results=results))
        
            
            
            """ For Starting and ending positions"""
        else:
            results = query_search.query.filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
            
            """ If nothing found redirect"""
            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))


            """ Run the function that parses the data from query object to python list object and assign them to sessions"""
            pop_data(results,variable)

            """ Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation"""
            flash("Please read the documentation for appropriate parameters", 'info')
            return redirect(url_for('results', title='Results'))
   
    else:

        """For Multi rsid or gene list, gene limited to just 3. if conditions not met
            redirects to search with a flashed message asking the user to check the format
        """
        if ',' in variable:
            try:  
                """ Creates a list of elements from user input seperated by comma : e.g:['a','b','c']"""
                variable = [i.strip() for i in variable.split(',')]

                rs_lst = []
                gene_lst=[]
                """ Assign the results to rs list if it has rs id format or other items from user input to gene list"""
                for i in variable:
                    """ Check rs id format"""
                    if i.startswith('rs') and isinstance(int(i[2:]),int):
                        rs_lst.append(i)
                    else:
                        """assigns everything else to gene list"""
                        gene_lst.append(i)
                

                """ Sanity check to make sure both lists have items concurrently or no values at the same time"""
                if len(rs_lst)>0 and len(gene_lst)>0:
                    raise Exception
                elif len(rs_lst)==0 and len(gene_lst)==0:
                    raise Exception
                else:
                    pass
                
                """ Not allowing for more than 3 genes"""
                if len(gene_lst)>3:
                    flash("Sorry Only 3 genes allowed", 'danger')
                    return redirect(url_for('search'))

                
                else:
                    """ If all checks passed and has rs_id, do query"""
                    if len(rs_lst) != 0:
                        results = query_search.query.filter(query_search.rs_val.in_(rs_lst)).all()

                        """ If no results found redirects to search to search something else"""
                        if not results:
                            flash("No result found, please search for another ID", 'info')
                            return redirect(url_for('search'))

                        
                        pop_data(results,rs_lst)
                        """ Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation"""
                        flash("Please read the documentation for appropriate parameters", 'info')
                        return redirect(url_for('results', title='Results'))
                        
                    else:
                        """ Automatically searches for the gene list if above condition is not met"""
                        results = query_search.query.filter(query_search.gene_name.in_(gene_lst)).all()
                       
                        """ If no results found redirects to search to search something else"""
                        if not results:
                            flash("No result found, please search for another ID", 'info')
                            return redirect(url_for('search'))
                            
                        """ Run the function that parses the data from query object to python list object and assign them to sessions"""
                        pop_data(results,gene_lst)

                        """ Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation"""
                        flash("Please read the documentation for appropriate parameters", 'info')
                        return redirect(url_for('results', title='Results')) 


            except Exception:
                """ If any of the above fails, assuming a bad user input"""
                flash("Sorry please check your format and try again", 'danger')
                return redirect(url_for('search'))


        
        elif variable.startswith('rs'):
            """ If its a single rs id value making sure no special characters are present (e.g. , . !)"""
            if variable.isalnum():
                results = query_search.query.filter(query_search.rs_val.like(variable)).all() 

                """ If no results found redirects to search to search something else"""
                if not results:
                    flash("No result found, please search for another ID", 'info')
                    return redirect(url_for('search'))

                """ Run the function that parses the data from query object to python list object and assign them to sessions"""
                pop_data(results,variable)

                """ Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation"""
                flash("Please read the documentation for appropriate parameters", 'info')
                return redirect(url_for('results', title='Results'))
            else:
                """ Bad characters in the input"""
                flash("Please Check your format", 'info')
                return redirect(url_for('search'))
        


        else:

            """ Only expecting a single gene"""
            variable = variable.upper()
            results = query_search.query.filter(query_search.gene_name.like(variable)).all()

            """ If no results found redirects to search to search something else"""
            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('search'))

            """ Run the function that parses the data from query object to python list object and assign them to sessions"""
            pop_data(results,variable)

            """ Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation"""
            flash("Please read the documentation for appropriate parameters", 'info')
            return redirect(url_for('results', title='Results', Results=results))



def convert_freq(pop):
    for i in pop:
        gf = ast.literal_eval(i['geno_freq'])
        gf_sum = sum(gf.values())
        gf_var = f"Hom-Ref:{round(gf['hom_ref']/gf_sum,2)}         Het:{round(gf['het']/gf_sum,2)}        Hom-Alt:{round(gf['hom_alt']/gf_sum,2)}"


        af = ast.literal_eval(i['allele_freq'])
        af_sum = sum(af.values())
        af_var = f"REF:{round(af['ref']/af_sum,2)}      ALT:{round(af['alt']/af_sum,2)}"

        i['geno_freq'] = gf_var
        i['allele_freq'] = af_var 

    return None





@app.route("/results", methods=['GET', 'POST'])
def results():
    try:
        results = json.loads(session['results'])
        gbr = json.loads(session['gbr'])
        jpt = json.loads(session['jpt'])
        mxl = json.loads(session['mxl'])
        pjl = json.loads(session['pjl'])
        yri = json.loads(session['yri'])
    except Exception:
        flash ('Please Search for SNPs first', 'info')
        return redirect(url_for('search'))


    convert_freq(gbr)
    convert_freq(jpt)
    convert_freq(mxl)
    convert_freq(pjl)
    convert_freq(yri)

    form = PopulationStatistics()
    if request.method == "POST":
        if form.validate_on_submit():
            if len(results) <= 1:
                flash("Can only perform statistics on two or more SNPs, please search for multiple SNPs", 'warning')
                return redirect(url_for('search'))
            return redirect(url_for('stats',pops=form.populations.data, stats=form.stats.data,bin=form.bin_size.data,step=form.step_size.data))
    return render_template('results.html', Results=results, GBR=gbr, JPT=jpt,MXL=mxl,PJL=pjl,YRI=yri, form=form)














""" Data for geneotyype array is stored as compressed data using similar technique to Run-length encoding (RLE)
    the compressed data includes list of lists [['2','a']]. This corresponds to [[0,0],[0,0]] : two samples that are
    homozygous alternative. The positions of the sample are also saved.  To decompress we use the value that corespond to 
    the compressed value and multiply with the integer
"""

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
@app.route("/stats/<pops>/<stats>/<bin>")
@app.route("/stats/<pops>/<stats>/<bin>/<step>")
def stats(pops= None, stats=None, bin=None, step = None):

    """If the data was store in sessions on a previous stats calc, load sessions so user can navigate around the site without having to redo the stats else, assumed to not have done any stats therefore expecting inputs for stats, population and bin size. If else also fails redirects results page to choose the parameters. 
    """
    if not pops:
        try:
            html_first_col=json.loads(session['overall_location'])
            html_gene=json.loads(session['gene_names'])
            all_stats=json.loads(session['all_stats'])
            pop_stats = json.loads(session['pop_stats'])
            all_fstat =json.loads(session['all_fstat'])
            nuc_div_plot1=json.loads(session['nuc_div_plot1'])
            nuc_div_plot2=json.loads(session['nuc_div_plot2'])
            hap_div_plot1=json.loads(session['hap_div_plot1'])
            hap_div_plot2=json.loads(session['hap_div_plot2'])
            taj_d_plot1=json.loads(session['taj_d_plot1'])
            taj_d_plot2=json.loads(session['taj_d_plot2'])
            gbr_fst_plt= json.loads(session['gbr_fst_plt'])
            jpt_fst_plt=json.loads(session['jpt_fst_plt'])
            mxl_fst_plt=json.loads(session['mxl_fst_plt'])
            pjl_fst_plt=json.loads(session['pjl_fst_plt'])
            yri_fst_plt=json.loads(session['yri_fst_plt'])
            
            return render_template('stats.html',
            html_first_col=html_first_col,
            html_gene=html_gene,
            all_stats=all_stats,
            pop_stats = pop_stats,
            all_fstat = all_fstat,
            nuc_div_plot1=nuc_div_plot1,
            nuc_div_plot2=nuc_div_plot2,
            hap_div_plot1=hap_div_plot1,
            hap_div_plot2=hap_div_plot2,
            taj_d_plot1=taj_d_plot1,
            taj_d_plot2=taj_d_plot2,
            gbr_fst_plt= gbr_fst_plt,
            jpt_fst_plt=jpt_fst_plt,
            mxl_fst_plt=mxl_fst_plt,
            pjl_fst_plt=pjl_fst_plt,
            yri_fst_plt=yri_fst_plt
            )
        except:
            pass


    else:

        """Expecting User to pass in valid values via GET request. Formatting issues can cause it to raise an error 
            and redirect them back to results page to properly pass them in. The user input format has to be specific 
            so they cant manually type them in the URL.
        """

        try:
            """ Checking if the user selected more than one population which is passed a string list i.e ("['pop1','pop2',]")"""
            if pops.startswith('[') :
                """ converts a string list into an acutal list type data
                str_lst = "['pop1','pop2']"
                type(str_lst) = str
                
                
                str_lst = "['pop1','pop2']"
                str_lst = ast.literal_eval(str_lst)
                type(str_lst) = lst
                """
                sel_pops = ast.literal_eval(pops)       
            else:
                """ Creating a list with one element if user selected just one population"""
                sel_pops = [pops]  
            
            """ Similarly checking if the user selected more than one stat which is passed a string list i.e ("['stat1','stat2',]")"""
            if stats.startswith('['):
                stats_sel = ast.literal_eval(stats) 
            else:
                stats_sel = [stats]

            """Get the bin size"""
            bin_size = int(bin)

            """Get the stepsize if user provided them"""
            if step:
                step_size = int(step)
            else:
                step_size =None
        except Exception:
            flash ('Please select the Stats and populations from this page', 'info')
            return redirect(url_for('results'))


        """Load the query session to perform calculations"""
        
        results = json.loads(session['results'])
        gen_pos = [int(i['pos']) for i in results]


        """Summary Stats for each population"""
        """ GBR"""
        if 'GBR' in sel_pops:
            """Load the data via sessions"""
            gbr = json.loads(session['gbr'])
            gbr_gt_data, gbr_freq = decompress(gbr)

            """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""        
            gbr_homo, gbr_nuc_div, gbr_hap_div, gbr_taj_d = gstat.get_main_stats(pop=gbr_gt_data, freq_data =gbr_freq, pos=gen_pos, stats=stats_sel)
            gbr_stats = ['GBR',gbr_homo, gbr_nuc_div, gbr_hap_div, gbr_taj_d]

            """Windowed Satats"""
            """PI"""
            gbr_win_pi= gstat.win_nuc_div(positions=gen_pos, pop=gbr_gt_data, bin_size=bin_size, step_size=step_size)

            """Tajimas D"""
            gbr_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=gbr_gt_data,bin_size=bin_size,step_size=step_size)

            """ Haplotype"""
            gbr_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=gbr_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            gbr = None



        """JPT"""
        """Load the data via sessions"""
        if 'JPT' in sel_pops:
            jpt = json.loads(session['jpt'])
            jpt_gt_data, jpt_freq = decompress(jpt)

            """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
            jpt_homo, jpt_nuc_div, jpt_hap_div, jpt_taj_d = gstat.get_main_stats(pop=jpt_gt_data, freq_data=jpt_freq, pos=gen_pos, stats=stats_sel)
            jpt_stats = ['JPT',jpt_homo, jpt_nuc_div, jpt_hap_div, jpt_taj_d]

            """Windowed Satats"""

            """PI"""
            jpt_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=jpt_gt_data,bin_size=bin_size,step_size=step_size)

            """Tajimas D"""
            jpt_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=jpt_gt_data,bin_size=bin_size,step_size=step_size)

            """ Haplotype"""
            jpt_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=jpt_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            jpt = None

        
        
        """MXL"""
        """Load the data via sessions"""
        if 'MXL' in sel_pops:
            mxl = json.loads(session['mxl'])
            mxl_gt_data, mxl_freq = decompress(mxl)
            """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
            mxl_homo, mxl_nuc_div, mxl_hap_div, mxl_taj_d = gstat.get_main_stats(pop=mxl_gt_data,freq_data =mxl_freq,pos=gen_pos,stats=stats_sel)
            mxl_stats = ['MXL',mxl_homo, mxl_nuc_div, mxl_hap_div, mxl_taj_d]

            """Windowed Satats"""
            """PI"""
            mxl_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=mxl_gt_data,bin_size=bin_size,step_size=step_size)

            """Tajimas D"""
            mxl_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=mxl_gt_data,bin_size=bin_size,step_size=step_size)

            """ Haplotype"""
            mxl_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=mxl_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            mxl = None



        """PJL"""
        """Load the data via sessions"""
        if 'PJL' in sel_pops:
            pjl = json.loads(session['pjl'])
            pjl_gt_data, pjl_freq = decompress(pjl)
            """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
            pjl_homo, pjl_nuc_div, pjl_hap_div, pjl_taj_d = gstat.get_main_stats(pop=pjl_gt_data,freq_data =pjl_freq,pos=gen_pos,stats=stats_sel)
            pjl_stats = ['PJL',pjl_homo, pjl_nuc_div, pjl_hap_div, pjl_taj_d]

            """Windowed Satats"""
            """PI"""
            pjl_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=pjl_gt_data,bin_size=bin_size,step_size=step_size)

            """Tajimas D"""
            pjl_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=pjl_gt_data,bin_size=bin_size,step_size=step_size)

            """ Haplotype"""
            pjl_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=pjl_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            pjl = None


        """YRI"""
        """Load the data via sessions"""
        if 'YRI' in sel_pops:
            yri = json.loads(session['yri'])
            yri_gt_data, yri_freq = decompress(yri)
            """Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D"""
            yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d = gstat.get_main_stats(pop=yri_gt_data,freq_data =yri_freq,pos=gen_pos,stats=stats_sel)
            yri_stats = ['YRI',yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d]

            """Windowed Satats"""
            """PI"""
            yri_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=yri_gt_data,bin_size=bin_size,step_size=step_size)

            """Tajimas D"""
            yri_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=yri_gt_data,bin_size=bin_size,step_size=step_size)

            """ Haplotype"""
            yri_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=yri_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            yri = None





        """For calculating fstats, creating a dictionary with key as the population code  and value as the array"""
        gt_dict = {}
        gt_freq = {}


        """This is for creating a dictionary of all the summary statistics calculated from above (no windows included)
            so its easier to use jija2 to dynamically create the tables and their stts for individual populations
        """
        pop_stats = {}


        """For creating plots"""
        plot_pi ={}
        plot_taj_d = {}
        plot_hap = {}

        """
        If the user seleceted the population assigning the genotype and frequency data to a diction with their key
        """
        if gbr:
            """For fst stats"""
            gt_dict['GBR'] = gbr_gt_data

            """Purely for calculating overall stats"""
            gt_freq['GBR'] = gbr_freq

            """For displauing on the table as individual stats"""
            pop_stats['GBR'] = gbr_stats

            """ For all plots"""
            plot_pi['GBR'] = gbr_win_pi
            plot_taj_d['GBR']=gbr_win_taj_D
            plot_hap['GBR'] = gbr_win_hap

        if jpt:
            """For fst stats"""
            gt_dict['JPT'] = jpt_gt_data

            """Purely for calculating overall stats"""
            gt_freq['JPT'] = jpt_freq

            """For displauing on the table as individual stats"""
            pop_stats['JPT'] = jpt_stats

            """ For all plots"""
            plot_pi['JPT'] = jpt_win_pi
            plot_taj_d['JPT']=jpt_win_taj_D
            plot_hap['JPT'] = jpt_win_hap

        if mxl:
            """For fst stats"""
            gt_dict['MXL'] = mxl_gt_data

            """Purely for calculating overall stats"""
            gt_freq['MXL'] = mxl_freq

            """For displauing on the table as individual stats"""
            pop_stats['MXL'] = mxl_stats

            """ For all plots"""
            plot_pi['MXL'] = mxl_win_pi
            plot_taj_d['MXL']=mxl_win_taj_D
            plot_hap['MXL'] = mxl_win_hap
        
        if pjl:
            """For fst stats"""
            gt_dict['PJL'] =pjl_gt_data

            """Purely for calculating overall stats"""
            gt_freq['PJL'] = pjl_freq

            """For displauing on the table as individual stats"""
            pop_stats['PJL'] = pjl_stats

            """ For all plots"""
            plot_pi['PJL'] = pjl_win_pi
            plot_taj_d['PJL']=pjl_win_taj_D
            plot_hap['PJL'] = pjl_win_hap


        if yri:
            """For fst stats"""
            gt_dict['YRI'] = yri_gt_data

            """Purely for calculating overall stats"""
            gt_freq['YRI'] = yri_freq

            """For displauing on the table as individual stats"""
            pop_stats['YRI'] = yri_stats

            """ For all plots"""
            plot_pi['YRI'] = yri_win_pi
            plot_taj_d['YRI']=yri_win_taj_D
            plot_hap['YRI'] = yri_win_hap
            

        """ Create Fst only if more than one population is selected and FST in stats using the  gt_dict dictionary"""
        if len(pops) >1 and 'FST' in stats_sel:
            all_fstat = gstat.get_fstat(paris=sel_pops, gt_dict = gt_dict)

            """ returns a list of tupes and converting them to a diction for easy of access downstream"""
            all_fstat = dict(all_fstat)
        else:
            all_fstat = None


        
        
        
        
        """ Get the overall location and the associated genes"""
        
        first_col = results[0]
        last_col = results[-1]


        html_first_col = f"CHR:{first_col['chrom']} S:{first_col['pos']}- E:{last_col['pos']}"
        html_gene = set([i['gene_name'] for i in results if i['gene_name'] != None])
        html_gene = ', '.join(html_gene)

        
        
        """ Overall stats by combining all sleected populations genotype data into a single 3d array
            using the gt_dict and gt_freq created for doing fst calculation to combine all the individual populations
            into a single superpopulation and calculating an overall summary stat
            gtd = genotype data
            cts = genotype counts
        """
        all_pops_gtd =[gt_dict[i] for i in sel_pops]
        all_pops_cts = [gt_freq[i] for i in sel_pops]

        all_pops_gtd = gstat.overall_stats_gtd(all_pops_gtd)
        all_pops_cts =gstat.overall_stats_cts(all_pops_cts)
        all_homo,all_nuc_div,all_hap_div,all_taj_d =gstat.get_main_stats(pop=all_pops_gtd,freq_data=all_pops_cts,pos=gen_pos,stats = stats_sel)

        """ Creating a dictionary so its easier to display in HTMl"""
        all_stats ={'Observed Homozygosity':all_homo,'Nucleotide Diversity(pi)':all_nuc_div,'Haplotide Diversity':all_hap_div,'Tajima D':all_taj_d}




        """ Get the population stats in a specific format for html they can be displayed separately
            this is a list of lists with each list starting with their population code and stats after the
            population code element
        """
        pop_stats = [pop_stats[i] for i in sel_pops]



        
        """Get the avg window size for plotting x axis by creating windows based on the window size and step
            size and then getting their avg.
        """

        x_axis = gstat.avg_win(gen_pos, size=bin_size,step=step_size)


        """Create plots based on the user selected stats"""
        if 'Nucleotide Diversity' in stats_sel:
            nuc_div_plot1,nuc_div_plot2 = gstat.plot_nuc_div(plot_pi,x_axis,sel_pops)
        else:
            nuc_div_plot1,nuc_div_plot2 = (None,None)
        
        if 'Haplotype Diversity' in stats_sel:
            hap_div_plot1,hap_div_plot2 = gstat.plot_win_hap(plot_hap,x_axis,sel_pops)
        else:
            hap_div_plot1,hap_div_plot2 = (None,None)
    
        if 'Tajimas D' in stats_sel:
            taj_d_plot1,taj_d_plot2 = gstat.plot_win_taj_d(plot_taj_d,x_axis,sel_pops)
        else:
            taj_d_plot1,taj_d_plot2 = (None,None)




        """For creating windowed fst plot we first get the overall fst for all combinations in a dictionary format and then split them into their respective population dictionaries, with each dictionary having a different pariwase comparisions as the key and their values being an array of fst in each window 
        """
        gbr_win_fst = {}
        jpt_win_fst = {}
        mxl_win_fst = {}
        pjl_win_fst = {}
        yri_win_fst = {}

        if len(pops) >1 and 'FST' in stats_sel:
            all_win_fstat = gstat.get_win_fstat(paris=sel_pops, gt_dict = gt_dict,pos=gen_pos,bin_size=bin_size,step_size=step_size)
            all_win_fstat = dict(all_win_fstat)
            for key,value in all_win_fstat.items():
                if key.startswith('GBR'):
                    gbr_win_fst[key] = value
                elif key.startswith('JPT'):
                    jpt_win_fst[key] = value
                elif key.startswith('MXL'):
                    mxl_win_fst[key] = value
                elif key.startswith('PJL'):
                    pjl_win_fst[key] = value
                elif key.startswith('YRI'):
                    yri_win_fst[key] = value
                else:
                    pass


        """ This will only create the plot if the above dictionaries are filled, else None is assigned to dynamically plot in the HTML jinja2 engine allowing for user selection plots."""
            
        if gbr_win_fst:
            gbr_fst_plt = gstat.plot_win_FST(pop_FST=gbr_win_fst,position=x_axis)
        else:
            gbr_fst_plt =None

        if jpt_win_fst:
            jpt_fst_plt = gstat.plot_win_FST(pop_FST=jpt_win_fst,position=x_axis)
        else:
            jpt_fst_plt =None

        if mxl_win_fst:
            mxl_fst_plt = gstat.plot_win_FST(pop_FST=mxl_win_fst,position=x_axis)
        else:
            mxl_fst_plt =None


        if pjl_win_fst:
            pjl_fst_plt = gstat.plot_win_FST(pop_FST=pjl_win_fst,position=x_axis)
        else:
            pjl_fst_plt =None


        if yri_win_fst:
            yri_fst_plt = gstat.plot_win_FST(pop_FST=yri_win_fst,position=x_axis)
        else:
            yri_fst_plt =None



        """Load everything into session to be used in the download function. If sessions are already populated replaces them.
        """
        session['overall_location'] = json.dumps(html_first_col)
        session['gene_names'] = json.dumps(html_gene)
        session['all_stats'] = json.dumps(all_stats)
        session['pop_stats'] = json.dumps(pop_stats)
        session['all_fstat'] = json.dumps(all_fstat)
        session['nuc_div_plot1'] = json.dumps(nuc_div_plot1)
        session['nuc_div_plot2'] = json.dumps(nuc_div_plot2)
        session['hap_div_plot1'] = json.dumps(hap_div_plot1)
        session['hap_div_plot2'] = json.dumps(hap_div_plot2)
        session['taj_d_plot1'] = json.dumps(taj_d_plot1)
        session['taj_d_plot2'] = json.dumps(taj_d_plot2)
        session['gbr_fst_plt'] = json.dumps(gbr_fst_plt)
        session['jpt_fst_plt'] = json.dumps(jpt_fst_plt)
        session['mxl_fst_plt'] = json.dumps(mxl_fst_plt)
        session['pjl_fst_plt'] = json.dumps(pjl_fst_plt)
        session['yri_fst_plt'] = json.dumps(yri_fst_plt)

        
        




        return render_template('stats.html',
        html_first_col=html_first_col,
        html_gene=html_gene,
        all_stats=all_stats,
        pop_stats = pop_stats,
        all_fstat = all_fstat,
        nuc_div_plot1=nuc_div_plot1,
        nuc_div_plot2=nuc_div_plot2,
        hap_div_plot1=hap_div_plot1,
        hap_div_plot2=hap_div_plot2,
        taj_d_plot1=taj_d_plot1,
        taj_d_plot2=taj_d_plot2,
        gbr_fst_plt= gbr_fst_plt,
        jpt_fst_plt=jpt_fst_plt,
        mxl_fst_plt=mxl_fst_plt,
        pjl_fst_plt=pjl_fst_plt,
        yri_fst_plt=yri_fst_plt
        )



    """ If user has did not have any previous stats calculations or did not pass in any values"""
    flash ('Select Stats and populations from this page first', 'info')
    return redirect(url_for('results'))









@app.route('/download')
def download():
    si = StringIO()
    fields = [
        'chrom',
        'rs_val',
        'pos',
        'gene_name',
        'ref_allele',
        'alt_allele',
        'GBR_geno_freq',
        'GBR_allele_freq',
        'MXL_geno_freq',
        'MXL_allele_freq',
        'JPT_geno_freq',
        'JPT_allele_freq',
        'PJL_geno_freq',
        'PJL_allele_freq',
        'YRI_geno_freq',
        'YRI_allele_freq'
        ]
    results = json.loads(session['results'])
    gbr = json.loads(session['gbr'])
    jpt = json.loads(session['jpt'])
    mxl = json.loads(session['mxl'])
    pjl = json.loads(session['pjl'])
    yri = json.loads(session['yri'])

    # Recreate the list of dictionaries, adding inside each dictionary the population key/value pairs
    # The final list of dicts inside the file will be "results_print", 
    # while we use "row_enriched" to enrich each dict
    results_print = []
    row_enriched = {}

    # Take each dict from the list
    for row in results:
        #Copy the dict into a new dict
        row_enriched = row.copy()
        
        # Loop inside each population dict to find the frequencies for this rs_id
        # when found, add in the enriched_row_dict
        for row_pop in gbr:
           if row["rs_val"] == row_pop["rs_val_id"]:
                row_enriched["GBR_geno_freq"] = row_pop["geno_freq"]
                row_enriched["GBR_allele_freq"] = row_pop["allele_freq"]

        for row_pop in mxl:
           if row["rs_val"] == row_pop["rs_val_id"]:
                row_enriched["MXL_geno_freq"] = row_pop["geno_freq"]
                row_enriched["MXL_allele_freq"] = row_pop["allele_freq"]
        
        for row_pop in jpt:
           if row["rs_val"] == row_pop["rs_val_id"]:
                row_enriched["JPT_geno_freq"] = row_pop["geno_freq"]
                row_enriched["JPT_allele_freq"] = row_pop["allele_freq"]

        for row_pop in pjl:
           if row["rs_val"] == row_pop["rs_val_id"]:
                row_enriched["PJL_geno_freq"] = row_pop["geno_freq"]
                row_enriched["PJL_allele_freq"] = row_pop["allele_freq"]

        for row_pop in yri:
           if row["rs_val"] == row_pop["rs_val_id"]:
                row_enriched["YRI_geno_freq"] = row_pop["geno_freq"]
                row_enriched["YRI_allele_freq"] = row_pop["allele_freq"]

        results_print.append(row_enriched)
    


    cw = csv.DictWriter(si, fieldnames=fields)
    cw.writeheader()
    for stats in results_print:
        cw.writerow(stats)
    output = make_response(si.getvalue())
    output.headers["Content-Disposition"] = "attachment; filename=query_results.csv"
    output.headers["Content-type"] = "text/csv"
    return output



@app.route('/download_stats')
def download_stats():
    html_first_col=json.loads(session['overall_location'])
    html_gene=json.loads(session['gene_names'])
    if not html_gene:
        html_gene = "None"
    all_stats=json.loads(session['all_stats'])

    pop_stats = json.loads(session['pop_stats'])
    all_fstats = json.loads(session['all_fstat'])
    # print(type(all_fstats))

    """ Combining Fstats with populations stats for the download format"""
    for key,value in all_fstats.items():
        """Checks if the fstat is a GBR comparision with other populations"""
        if key.startswith('GBR'):
            temp = f"{key}:{value}"
            for i in pop_stats:
                """ Checks if the pop_stats element belonds to GBR and appends the string"""
                if i[0] == 'GBR':
                    i.append(temp)
                    break

            """Checks if the fstat is a JPT comparision with other populations"""
        elif key.startswith('JPT'):
            temp = f"{key}:{value}"
            for i in pop_stats:
                """ Checks if the pop_stats element belonds to JPT and appends the string"""
                if i[0] == 'JPT':
                    i.append(temp)
                    break
        elif key.startswith('MXL'):
            temp = [f"{key}:{value}"]
            for i in pop_stats:
                if i[0] == 'MXL':
                    i.extend(temp)
                    break
        elif key.startswith('PJL'):
            temp = f"{key}:{value}"
            for i in pop_stats:
                if i[0] == 'PJL':
                    i.append(temp)
                    break
        elif key.startswith('YRI'):
            temp = f"{key}:{value}"
            for i in pop_stats:
                if i[0] == 'YRI':
                    i.append(temp)
                    break
        else:
            pass

    """ Create a string IO object"""
    si = StringIO()
    cw = csv.writer(si,delimiter='\t')

    """ String io needs the values to be in a list"""
    first_col_header = ["LOCATION","GENES","Homzygosity","Nucleotide Diversity","Haplotype Diversity","Tajimas D"]
    cw.writerow(first_col_header)

    """ From all_stats each item is store as a list and indexes correspond to the stats retrieval order as follows:
        Homzygosity, Nucleotide Diversity, Haplotype Diversity, Tajimas D
    """
    first_col_values = [html_first_col, html_gene,all_stats['Observed Homozygosity'],all_stats['Nucleotide Diversity(pi)'],all_stats['Haplotide Diversity'],all_stats['Tajima D']]
    cw.writerow(first_col_values)

    """ Create a space between the """
    cw.writerow('')
    cw.writerow('')
    cw.writerow('')

    pop_headers = ["Population","Homzygosity","Nucleotide Diversity","Haplotype Diversity","Tajimas D","FST"]
    cw.writerow(pop_headers)
    for i in pop_stats:
        cw.writerow(i)
    
    mem = BytesIO()
    mem.write(si.getvalue().encode())
    mem.seek(0)
    si.close()
    return send_file(
        mem,
        as_attachment=True,
        attachment_filename='test.txt',
        mimetype='text/plain'
    )








@app.route("/help", methods=['GET', 'POST'])
def help():
    return render_template('help.html', title='Contact')







# Create Custom Error Pages
# Invalid URL
@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404


# Internal Server Error
@app.errorhandler(500)
def page_not_found(e):
    return render_template("500.html"), 500
