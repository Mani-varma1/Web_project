from flask import Flask, render_template
import plotly.graph_objects as go
import plotly
import json
import plotly.express as px

app = Flask(__name__)


@app.route('/')
def index():
    # here is the data from the stats
    # the data needs to be converted into a dictionary for the selected populations
    # the positions needs to have an average

    # eg Tajimas D
    TD = {'GBR': [0., 0., 0., 0.001, 0.], 'JPT': [0., 0., 0., 0., 0.0, ], 'MXL': [0., 0., 0.01, 0.0, 0.5],
          'PJL': [0.0, 0, 0.002, 0.0, 0.005], 'YRI': [0.0, 0.0, 0., 0.1, 0.0]}

    # eg Halploid diversity
    HD = {'GBR': [0., 0., 0., 0.001, 0.1], 'JPT': [0., 0., 0., 0., 0.2], 'MXL': [0., 0., 0.01, 0.0, 0.5],
          'PJL': [0.0, 0, 0.002, 0.0, 0.005], 'YRI': [0.0, 0.0, 0., 0.1, 0.0]}

    # eg Nucleotide diversity
    ND = {'GBR': [0., 0., 0., 0.001, 0.1], 'JPT': [0., 0., 0., 0., 0.2], 'MXL': [0., 0., 0.01, 0.0, 0.3],
          'PJL': [0.0, 0, 0.002, 0.0, 0.4], 'YRI': [0.0, 0.0, 0., 0.1, 0.5]}

    # eg average positions
    position = [16050075, 16050375, 16050425, 16050825, 16051075]

    user_selected_stats = ['TAJ_D', 'Haplotide', 'Nucleotide']



#creating graphs

    if 'TAJ_D' not in user_selected_stats:
        graph1JSON = None
        graph1JSON_CHECK = None
        graph2JSON = None
        graph2JSON_CHECK = None

    else:
        #creating Tajima's D graph
        fig_1 = go.Figure()
        buttons = []
        i = 0

        #iterating through dictionary and addding each population to graph
        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_1.add_trace(obj)

            #args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(TD)
            args[i] = True

            #creating button object for each population
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

            #add button to list
            buttons.append(button)
            i += 1

        #adding buttons
        fig_1.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])
        #adding axis names
        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

        #graph containing all populations without buttons
        fig_2 = go.Figure()
        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)

        fig_2.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph1JSON_CHECK = 'check'
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON_CHECK = 'check'

    if 'Haplotide' not in user_selected_stats:
        graph3JSON = None
        graph3JSON_CHECK = None
        graph4JSON = None
        graph4JSON_CHECK = None

    else:
        #creating Haplotype graphs
        fig_3 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary and adding each population
        for x in HD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_3.add_trace(obj)

        #args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(HD)
            args[i] = True

        #create button object for each pop
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

        #add button to list
            buttons.append(button)
            i += 1

        #add buttons
        fig_3.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])
        #add axis names
        fig_3.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

    # graph containing all populations without buttons
        fig_4 = go.Figure()
        for x in HD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_4.add_trace(obj)

        fig_4.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

        graph3JSON = json.dumps(fig_3, cls=plotly.utils.PlotlyJSONEncoder)
        graph3JSON_CHECK = 'check'
        graph4JSON = json.dumps(fig_4, cls=plotly.utils.PlotlyJSONEncoder)
        graph4JSON_CHECK = 'check'


    if 'Nucleotide' not in user_selected_stats:
        graph5JSON = None
        graph5JSON_CHECK = None
        graph6JSON = None
        graph6JSON_CHECK = None

    else:
        # Creating nucleotide Diversity graph
        fig_5 = go.Figure()
        buttons = []
        i = 0

        #iterating through dictionary and adding each population
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_5.add_trace(obj)

        ##args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(ND)
            args[i] = True

        #creating button for each population
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

            # adding button
            buttons.append(button)
            i += 1

        fig_5.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        #adding axis names
        fig_5.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Nucleotide Diversity")

     #graph for all populations without buttons
        fig_6 = go.Figure()
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_6.add_trace(obj)
            fig_6.update_layout(
                xaxis_title='Position (Base pairs)',
                yaxis_title="Nucleotide Diversity")

        graph5JSON = json.dumps(fig_5, cls=plotly.utils.PlotlyJSONEncoder)
        graph5JSON_CHECK = 'check'
        graph6JSON = json.dumps(fig_6, cls=plotly.utils.PlotlyJSONEncoder)
        graph6JSON_CHECK = 'check'

    return render_template('plotly.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON,
                           graph3JSON=graph3JSON, graph4JSON=graph4JSON,
                           graph5JSON=graph5JSON, graph6JSON=graph6JSON, graph1JSON_CHECK = graph1JSON_CHECK,
                           graph2JSON_CHECK = graph2JSON_CHECK, graph3JSON_CHECK = graph3JSON_CHECK,
                           graph4JSON_CHECK = graph4JSON_CHECK, graph5JSON_CHECK = graph5JSON_CHECK,
                           graph6JSON_CHECK = graph6JSON_CHECK)

if __name__ == '__main__':
	app.run(debug=True)