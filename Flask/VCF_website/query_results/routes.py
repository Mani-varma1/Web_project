from flask import render_template, url_for, flash, redirect, request, session,Blueprint,make_response
from VCF_website.models import snp_GBR,snp_JPT,snp_MXL,snp_PJL,snp_YRI,query_search
from VCF_website.query_results.results_utils import pop_data, convert_freq
from VCF_website.query_results.forms import PopulationStatistics
from io import StringIO
import json
import csv

query_results = Blueprint("query_results",__name__)

@query_results.route("/loading", methods=['GET', 'POST'])
def loading(search):
    variable = search
    # Explicitly clear sessions
    session.clear()
    
    
    # Check if the user input was using Chromosome locations which passed to this route as a dictionary
    if isinstance(variable, dict):
        
        # For single location
        if variable["end_pos"] == None:
            # queries based on starting position and chr
            results = query_search.query.filter(query_search.pos.like(variable['start_pos'])).filter(query_search.chrom == '{}'.format(variable["chr"])).all()

            # If nothing found redirect back to search with a flashed message notifying the user that nothing was found
            if not results:
                flash("No results found, please search for another ID", 'info')
                return redirect(url_for('main.search'))

            # Run the function that parses the data from query object to python list object and assign them to sessions from results_utils.py
            pop_data(results,variable)

            # Flash message so user can have a better understanding of the appropriate window size and step sizes inputs by reading the documentation
            flash("Please read the documentation for appropriate parameters", 'info')            
            return redirect(url_for('query_results.results', title='Results'))
        
            
            
        else:
            if int(variable['end_pos']) >= 16050075:
                # For a range between Starting and ending positions base pair queries based on starting position and chr
                results = query_search.query.filter(query_search.pos >= int(variable['start_pos'])).filter(query_search.pos <= int(variable['end_pos'])).filter(query_search.chrom == '{}'.format(variable['chr'])).all()
                
                # If nothing found redirect and inform
                if not results:
                    flash("No result found, please search for another ID", 'info')
                    return redirect(url_for('main.search'))


                # Run the function that parses the data from query object to python list object and assign them to sessions
                pop_data(results,variable)

                # Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation
                flash("Please read the documentation for appropriate parameters", 'info')
                return redirect(url_for('query_results.results', title='Results'))
            else:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('main.search'))
   
    else:
        # Both rsid and gene form submission are recieved as a string of a single value e.g. "rs123" or "GENE1" or as a string with commas
        # e.g. "rs123,rs456,...." or "GENE1,GENE2,...". The formating is strict and if other unexpected characters provided, should redirect 
        # back to search url with flashed message informing user of the format issue.

        # For Multi rsid or gene list, gene limited to just 3. if conditions not met
        # redirects to search with a flashed message asking the user to check the format
        if ',' in variable:
            try:  
                # Creates a list of elements from user input seperated by comma : e.g:['a','b','c']
                variable = [i.strip() for i in variable.split(',')]

                rs_lst = []
                gene_lst=[]
                # Assign the results to rs list if it has rs id format or other items from user input to gene list
                for i in variable:
                    # Check rs id format
                    if i.startswith('rs') and isinstance(int(i[2:]),int):
                        rs_lst.append(i)
                    else:
                        # assigns everything else to gene list
                        gene_lst.append(i.upper().strip())
                

                # Sanity check to make sure both lists have items concurrently or no values at the same time
                if len(rs_lst)>0 and len(gene_lst)>0:
                    raise Exception
                elif len(rs_lst)==0 and len(gene_lst)==0:
                    raise Exception
                else:
                    pass
                
                # Not allowing for more than 3 genes
                if len(gene_lst)>3:
                    flash("Sorry Only 3 genes allowed", 'danger')
                    return redirect(url_for('main.search'))

                
                else:
                    # If all checks passed and has rs_id, do query
                    if len(rs_lst) != 0:
                        results = query_search.query.filter(query_search.rs_val.in_(rs_lst)).all()

                        # If no results found redirects to search to search something else
                        if not results:
                            flash("No result found, please search for another ID", 'info')
                            return redirect(url_for('main.search'))

                        
                        pop_data(results,rs_lst)
                        # Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation
                        flash("Please read the documentation for appropriate parameters", 'info')
                        return redirect(url_for('query_results.results', title='Results'))
                        
                    else:
                        # Automatically searches for the gene list if above condition is not met
                        results = query_search.query.filter(query_search.gene_name.in_(gene_lst)).all()
                       
                        # If no results found redirects to search to search something else
                        if not results:
                            flash("No result found, please search for another ID", 'info')
                            return redirect(url_for('main.search'))
                            
                        # Run the function that parses the data from query object to python list object and assign them to sessions
                        pop_data(results,gene_lst)

                        # Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation
                        flash("Please read the documentation for appropriate parameters", 'info')
                        return redirect(url_for('query_results.results', title='Results')) 


            except Exception:
                # If any of the above fails, assuming a bad user input
                flash("Sorry please check your format and try again", 'danger')
                return redirect(url_for('main.search'))

   
            # If comma not in the list, expecting all the values are single so if it starts with rs, checking for dbSNP search.
        elif variable.startswith('rs'):
            # If its a single rs id value making sure no special characters are present (e.g: $ . !)
            if variable.isalnum():
                results = query_search.query.filter(query_search.rs_val.like(variable)).all() 

                # If no results found redirects to search to search something else
                if not results:
                    flash("No result found, please search for another ID", 'info')
                    return redirect(url_for('main.search'))

                # Run the function that parses the data from query object to python list object and assign them to sessions
                pop_data(results,variable)

                # Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation
                flash("Please read the documentation for appropriate parameters", 'info')
                return redirect(url_for('query_results.results', title='Results'))
            else:
                # Bad characters in the input there asking for users to check the format
                flash("Please Check your format", 'info')
                return redirect(url_for('main.search'))
        


        else:
            # Some Genes have - in them to for multiple loci?

            # Only expecting a single gene
            variable = variable.upper()
            results = query_search.query.filter(query_search.gene_name.like(variable)).all()

            # If no results found redirects to search to search something else
            if not results:
                flash("No result found, please search for another ID", 'info')
                return redirect(url_for('main.search'))

            # Run the function that parses the data from query object to python list object and assign them to sessions
            pop_data(results,variable)

            # Flash message so user can have a better understanding of the input window size and step sizes by reading the documentation
            flash("Please read the documentation for appropriate parameters", 'info')
            return redirect(url_for('query_results.results', title='Results'))


@query_results.route("/results", methods=['GET', 'POST'])
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
        return redirect(url_for('main.search'))


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
                return redirect(url_for('main.search'))
            return redirect(url_for('statistics.stats',pops=form.populations.data, stats=form.stats.data,bin=form.bin_size.data,step=form.step_size.data))
    return render_template('results.html', Results=results, GBR=gbr, JPT=jpt,MXL=mxl,PJL=pjl,YRI=yri, form=form)


@query_results.route('/download')
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
