from flask import render_template, url_for, flash, redirect, request, session, make_response,send_file, Blueprint
import VCF_website.statistics.genome_stats as gstat
from VCF_website.statistics.stats_utils import decompress
import ast
import json
from io import StringIO,BytesIO
import csv
from werkzeug.wrappers import Response

statistics = Blueprint("statistics",__name__)

@statistics.route("/stats")
@statistics.route("/stats/<pops>/<stats>/<bin>")
@statistics.route("/stats/<pops>/<stats>/<bin>/<step>")
def stats(pops= None, stats=None, bin=None, step = None):

    # If the data was store in sessions on a previous stats calc, load sessions so user can navigate around the site 
    # without having to redo the stats else, assumed to not have done any stats therefore expecting inputs for stats, population and bin size. 
    # If else also fails redirects results page to choose the parameters. 
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

        # Expecting User to pass in valid values via GET request. Formatting issues can cause it to raise an error 
        # and redirect them back to results page to properly pass them in. The user input format has to be specific 
        # so they cant manually type them in the URL.
        try:
            # Checking if the user selected more than one population which is passed a string list i.e ("['pop1','pop2',]")
            if pops.startswith('[') :
                # converts a string list into an acutal list type data
                sel_pops = ast.literal_eval(pops)       
            else:
                """ Creating a list with one element if user selected just one population"""
                sel_pops = [pops]  
            
            # Similarly checking if the user selected more than one stat which is passed a string list i.e ("['stat1','stat2',]")
            if stats.startswith('['):
                stats_sel = ast.literal_eval(stats) 
            else:
                stats_sel = [stats]

            # Get the bin size
            bin_size = int(bin)

            # Get the stepsize if user provided them
            if step:
                step_size = int(step)
            else:
                step_size =None
        except Exception:
            flash ('Please select the Stats and populations from this page', 'info')
            return redirect(url_for('query_results.results'))


        # Load the query session to perform calculations    
        results = json.loads(session['results'])
        gen_pos = [int(i['pos']) for i in results]


        # Summary Stats for each population
        #  GBR
        if 'GBR' in sel_pops:
            # Load the data via sessions
            gbr = json.loads(session['gbr'])
            gbr_gt_data, gbr_freq = decompress(gbr)

            # Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D        
            gbr_homo, gbr_nuc_div, gbr_hap_div, gbr_taj_d = gstat.get_main_stats(pop=gbr_gt_data, freq_data =gbr_freq, pos=gen_pos, stats=stats_sel)
            gbr_stats = ['GBR',gbr_homo, gbr_nuc_div, gbr_hap_div, gbr_taj_d]

            # Windowed Satats
            # PI
            gbr_win_pi= gstat.win_nuc_div(positions=gen_pos, pop=gbr_gt_data, bin_size=bin_size, step_size=step_size)

            # Tajimas D
            gbr_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=gbr_gt_data,bin_size=bin_size,step_size=step_size)

            #  Haplotype
            gbr_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=gbr_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            gbr = None


        # JPT
        # Load the data via sessions
        if 'JPT' in sel_pops:
            jpt = json.loads(session['jpt'])
            jpt_gt_data, jpt_freq = decompress(jpt)

            # Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D
            jpt_homo, jpt_nuc_div, jpt_hap_div, jpt_taj_d = gstat.get_main_stats(pop=jpt_gt_data, freq_data=jpt_freq, pos=gen_pos, stats=stats_sel)
            jpt_stats = ['JPT',jpt_homo, jpt_nuc_div, jpt_hap_div, jpt_taj_d]

            # Windowed Satats

            # PI
            jpt_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=jpt_gt_data,bin_size=bin_size,step_size=step_size)

            # Tajimas D
            jpt_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=jpt_gt_data,bin_size=bin_size,step_size=step_size)

            # Haplotype
            jpt_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=jpt_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            jpt = None
        
        
        # MXL
        # Load the data via sessions
        if 'MXL' in sel_pops:
            mxl = json.loads(session['mxl'])
            mxl_gt_data, mxl_freq = decompress(mxl)
            # Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D
            mxl_homo, mxl_nuc_div, mxl_hap_div, mxl_taj_d = gstat.get_main_stats(pop=mxl_gt_data,freq_data =mxl_freq,pos=gen_pos,stats=stats_sel)
            mxl_stats = ['MXL',mxl_homo, mxl_nuc_div, mxl_hap_div, mxl_taj_d]

            # Windowed Satats
            # PI
            mxl_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=mxl_gt_data,bin_size=bin_size,step_size=step_size)

            # Tajimas D
            mxl_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=mxl_gt_data,bin_size=bin_size,step_size=step_size)

            # Haplotype
            mxl_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=mxl_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            mxl = None


        # PJL
        # Load the data via sessions
        if 'PJL' in sel_pops:
            pjl = json.loads(session['pjl'])
            pjl_gt_data, pjl_freq = decompress(pjl)
            # Get homozygositym nucleotide diversity, haplotype diversity, and Tajimas D
            pjl_homo, pjl_nuc_div, pjl_hap_div, pjl_taj_d = gstat.get_main_stats(pop=pjl_gt_data,freq_data =pjl_freq,pos=gen_pos,stats=stats_sel)
            pjl_stats = ['PJL',pjl_homo, pjl_nuc_div, pjl_hap_div, pjl_taj_d]

            # Windowed Satats
            # PI
            pjl_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=pjl_gt_data,bin_size=bin_size,step_size=step_size)

            # Tajimas D
            pjl_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=pjl_gt_data,bin_size=bin_size,step_size=step_size)

            # Haplotype
            pjl_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=pjl_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            pjl = None


        # YRI
        # Load the data via sessions
        if 'YRI' in sel_pops:
            yri = json.loads(session['yri'])
            yri_gt_data, yri_freq = decompress(yri)
            # Get homozygosity nucleotide diversity, haplotype diversity, and Tajimas D
            yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d = gstat.get_main_stats(pop=yri_gt_data,freq_data =yri_freq,pos=gen_pos,stats=stats_sel)
            yri_stats = ['YRI',yri_homo,yri_nuc_div,yri_hap_div,yri_taj_d]

            # Windowed Stats
            # PI
            yri_win_pi=  gstat.win_nuc_div(positions=gen_pos,pop=yri_gt_data,bin_size=bin_size,step_size=step_size)

            # Tajimas D
            yri_win_taj_D= gstat.win_tajima_d(positions=gen_pos,pop=yri_gt_data,bin_size=bin_size,step_size=step_size)

            # Haplotype
            yri_win_hap = gstat.win_haplotype_div(positions=gen_pos,pop=yri_gt_data,bin_size=bin_size,step_size=step_size)
        else:
            yri = None



        # For calculating fstats, creating a dictionary with key as the population code  and value as the array
        gt_dict = {}
        gt_freq = {}


        # This is for creating a dictionary of all the summary statistics calculated from above (no windows included)
        # so its easier to use jija2 to dynamically create the tables and their stts for individual populations
        pop_stats = {}


        #For creating plots
        plot_pi ={}
        plot_taj_d = {}
        plot_hap = {}

        
        # If the user seleceted the population assigning the genotype and frequency data to a diction with their key
        if gbr:
            # For fst stats
            gt_dict['GBR'] = gbr_gt_data

            # Purely for calculating overall stats
            gt_freq['GBR'] = gbr_freq

            # For displauing on the table as individual stats
            pop_stats['GBR'] = gbr_stats

            # For all plots
            plot_pi['GBR'] = gbr_win_pi
            plot_taj_d['GBR']=gbr_win_taj_D
            plot_hap['GBR'] = gbr_win_hap

        if jpt:
            # For fst stats
            gt_dict['JPT'] = jpt_gt_data

            # Purely for calculating overall stats
            gt_freq['JPT'] = jpt_freq

            # For displauing on the table as individual stats
            pop_stats['JPT'] = jpt_stats

            # For all plots
            plot_pi['JPT'] = jpt_win_pi
            plot_taj_d['JPT']=jpt_win_taj_D
            plot_hap['JPT'] = jpt_win_hap

        if mxl:
            # For fst stats
            gt_dict['MXL'] = mxl_gt_data

            # Purely for calculating overall stats
            gt_freq['MXL'] = mxl_freq

            # For displauing on the table as individual stats
            pop_stats['MXL'] = mxl_stats

            # For all plots
            plot_pi['MXL'] = mxl_win_pi
            plot_taj_d['MXL']=mxl_win_taj_D
            plot_hap['MXL'] = mxl_win_hap
        
        if pjl:
            # For fst stats
            gt_dict['PJL'] =pjl_gt_data

            # Purely for calculating overall stats
            gt_freq['PJL'] = pjl_freq

            # For displauing on the table as individual stats
            pop_stats['PJL'] = pjl_stats

            # For all plots
            plot_pi['PJL'] = pjl_win_pi
            plot_taj_d['PJL']=pjl_win_taj_D
            plot_hap['PJL'] = pjl_win_hap


        if yri:
            # For fst stats
            gt_dict['YRI'] = yri_gt_data

            # Purely for calculating overall stats
            gt_freq['YRI'] = yri_freq

            # For displauing on the table as individual stats
            pop_stats['YRI'] = yri_stats

            # For all plots
            plot_pi['YRI'] = yri_win_pi
            plot_taj_d['YRI']=yri_win_taj_D
            plot_hap['YRI'] = yri_win_hap
            

        # Create Fst only if more than one population is selected and FST in stats using the  gt_dict dictionary
        if len(pops) >1 and 'FST' in stats_sel:
            all_fstat = gstat.get_fstat(paris=sel_pops, gt_dict = gt_dict)

            # returns a list of tupes and converting them to a diction for easy of access downstream
            all_fstat = dict(all_fstat)
        else:
            all_fstat = None


        
        # Get the overall location and the associated genes
        first_col = results[0]
        last_col = results[-1]

        html_first_col = f"CHR:{first_col['chrom']} S:{first_col['pos']}- E:{last_col['pos']}"
        html_gene = set([i['gene_name'] for i in results if i['gene_name'] != None])
        html_gene = ', '.join(html_gene)
        
        
        # Overall stats by combining all sleected populations genotype data into a single 3d array
        # using the gt_dict and gt_freq created for doing fst calculation to combine all the individual populations
        # into a single superpopulation and calculating an overall summary stat
        # gtd = genotype data
        # cts = genotype counts
        all_pops_gtd =[gt_dict[i] for i in sel_pops]
        all_pops_cts = [gt_freq[i] for i in sel_pops]

        all_pops_gtd = gstat.overall_stats_gtd(all_pops_gtd)
        all_pops_cts =gstat.overall_stats_cts(all_pops_cts)
        all_homo,all_nuc_div,all_hap_div,all_taj_d =gstat.get_main_stats(pop=all_pops_gtd,freq_data=all_pops_cts,pos=gen_pos,stats = stats_sel)


        # Creating a dictionary so its easier to display in HTMl"""
        all_stats ={'Observed Homozygosity':all_homo,'Nucleotide Diversity(pi)':all_nuc_div,'Haplotide Diversity':all_hap_div,'Tajima D':all_taj_d}


        # Get the population stats in a specific format for html they can be displayed separately
        # this is a list of lists with each list starting with their population code and stats after the
        # population code element
        pop_stats = [pop_stats[i] for i in sel_pops]


        # Get the avg window size for plotting x axis by creating windows based on the window size and step
        # size and then getting their avg.
        x_axis = gstat.avg_win(gen_pos, size=bin_size,step=step_size)


        # Create plots based on the user selected stats
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


        # For creating windowed fst plot we first get the overall fst for all combinations in a dictionary format 
        # and then split them into their respective population dictionaries, with each dictionary having a different 
        # pariwase comparisions as the key and their values being an array of fst in each window 
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


        # This will only create the plot if the above dictionaries are filled, 
        # else None is assigned to dynamically plot in the HTML jinja2 
        # engine allowing for user selection plots.
            
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


        # Load everything into session to be used in the download function. 
        # If sessions are already populated replaces them.
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


    # If user has did not have any previous stats calculations or did not pass in any values
    flash ('Select Stats and populations from this page first', 'info')
    return redirect(url_for('query_results.results'))



@statistics.route('/download_stats')
def download_stats():
    html_first_col=json.loads(session['overall_location'])
    html_gene=json.loads(session['gene_names'])
    if not html_gene:
        html_gene = "None"
    all_stats=json.loads(session['all_stats'])

    pop_stats = json.loads(session['pop_stats'])
    all_fstats = json.loads(session['all_fstat'])
    # print(type(all_fstats))

    # Combining Fstats with populations stats for the download format
    if all_fstats:
        for key,value in all_fstats.items():
            # Checks if the fstat is a GBR comparision with other populations
            if key.startswith('GBR'):
                temp = f"{key}:{value}"
                for i in pop_stats:
                    # Checks if the pop_stats element belonds to GBR and appends the string
                    if i[0] == 'GBR':
                        i.append(temp)
                        break

                # Checks if the fstat is a JPT comparision with other populations
            elif key.startswith('JPT'):
                temp = f"{key}:{value}"
                for i in pop_stats:
                    # Checks if the pop_stats element belonds to JPT and appends the string
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
                
    # Create a string IO object
    si = StringIO()
    cw = csv.writer(si,delimiter='\t')

    # String io needs the values to be in a list
    first_col_header = ["LOCATION","GENES","Homzygosity","Nucleotide Diversity","Haplotype Diversity","Tajimas D"]
    cw.writerow(first_col_header)

    # From all_stats each item is store as a list and indexes correspond to the stats retrieval order as follows:
    # Homzygosity, Nucleotide Diversity, Haplotype Diversity, Tajimas D
    first_col_values = [html_first_col, html_gene,all_stats['Observed Homozygosity'],all_stats['Nucleotide Diversity(pi)'],all_stats['Haplotide Diversity'],all_stats['Tajima D']]
    cw.writerow(first_col_values)

    # Create a space between the 
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
