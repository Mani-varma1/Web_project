from flask import Flask, render_template
import plotly.graph_objects as go
import plotly
import plotly.express as px
import json



app = Flask(__name__)


@app.route('/')
def index():

# here is the data from the stats
# the data needs to be converted into a dictionary for the selected populations
# the positions needs to have an average

# eg Tajimas D
    #TD = {}
    TD = {'GBR': [0., 0., 0., 0.001, 0.], 'JPT': [0., 0., 0., 0., 0.0, ], 'MXL': [0., 0., 0.01, 0.0, 0.5],
      'PJL': [0.0, 0, 0.002, 0.0, 0.005], 'YRI': [0.0, 0.0, 0., 0.1, 0.0]}

# eg Halploid diversity
    #HD = {}
    HD = {'GBR': [0., 0., 0.,0.001, 0.1], 'JPT': [0., 0., 0., 0., 0.2], 'MXL': [0., 0., 0.01, 0.0, 0.5],
      'PJL': [0.0, 0, 0.002, 0.0, 0.005], 'YRI': [0.0, 0.0, 0., 0.1, 0.0]}

#eg Nucleotide diversity
    #ND = {}
    ND = {'GBR': [0., 0., 0., 0.001, 0.1], 'JPT': [0., 0., 0., 0., 0.2], 'MXL': [0., 0., 0.01, 0.0, 0.3],
      'PJL': [0.0, 0, 0.002, 0.0, 0.4], 'YRI': [0.0, 0.0, 0., 0.1, 0.5]}

# average positions
    position = [16050075, 16050375, 16050425, 16050825, 16051075]

# graphs for all three stats
    if len(HD) != 0 and len(ND) != 0 and len(TD) != 0:

        #creating Tajimas D graph indivual
        fig_1 = go.Figure()
        buttons = []
        i = 0

        #iterating through dict and adding each pop to graph
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
        # adding axis names
        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

    # graph for all populations Tajima D without buttons
        fig_2 = go.Figure()

        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)

        fig_2.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

    # Creating Halploid Diversity graph individual

        fig_3 = go.Figure()
        buttons = []
        i = 0

    # iterating through dictionary
        for x in HD.items():
        # adding scatter line for each pop
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_3.add_trace(obj)

        #args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(HD)
            args[i] = True

        # create button object for each pop
            button = dict(label=x[0],
                        method="update",
                        args=[{"visible": args}])

        # add button to list
            buttons.append(button)
            i += 1

        fig_3.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_3.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

    # graph for overall Haploid Diversity_
        fig_4 = go.Figure()
        for x in HD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_4.add_trace(obj)

        fig_4.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

    # Creating nucleotide Diversity graph individual

        fig_5 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_5.add_trace(obj)
            args = [False] * len(ND)
            args[i] = True
            button = dict(label=x[0],
                        method="update",
                        args=[{"visible": args}])
            buttons.append(button)
            i += 1

        fig_5.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_5.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Nucleotide Diversity")

        # graph for overall _
        fig_6 = go.Figure()
        for x in HD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_6.add_trace(obj)
            fig_6.update_layout(
                xaxis_title='Position (Base pairs)',
                yaxis_title="Nucleotide Diversity")

        #coverting to json objects
        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)
        graph3JSON = json.dumps(fig_3, cls=plotly.utils.PlotlyJSONEncoder)
        graph4JSON = json.dumps(fig_4, cls=plotly.utils.PlotlyJSONEncoder)
        graph5JSON = json.dumps(fig_5, cls=plotly.utils.PlotlyJSONEncoder)
        graph6JSON = json.dumps(fig_6, cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('stats_1.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON,
                       graph3JSON=graph3JSON, graph4JSON=graph4JSON,
                       graph5JSON=graph5JSON, graph6JSON=graph6JSON)

   ########################################################################################## only Tajimas D and Haploid Diversity graphs##################################################################
    if len(TD) != 0 and len(HD) !=0:

        # creating Tajimas D graph individual
        fig_1 = go.Figure()
        buttons = []
        i = 0

        # iterating through dict
        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_1.add_trace(obj)

            # args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(TD)
            args[i] = True

            # creating button object for each population
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

            # add button to list
            buttons.append(button)
            i += 1

        fig_1.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])
        # adding axis names
        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

        # graph for all populations Tajima D
        fig_2 = go.Figure()

        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)

        fig_2.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

        # Creating Halploid Diversity graph individual

        fig_3 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary
        for x in HD.items():
            # adding scatter line for each pop
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_3.add_trace(obj)

            # args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(HD)
            args[i] = True

            # create button object for each pop
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

            # add button to list
            buttons.append(button)
            i += 1

        fig_3.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_3.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

        # graph for overall Haploid Diversity_
        fig_4 = go.Figure()
        for x in HD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_4.add_trace(obj)
        fig_4.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)
        graph3JSON = json.dumps(fig_3, cls=plotly.utils.PlotlyJSONEncoder)
        graph4JSON = json.dumps(fig_4, cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('stats_2.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON,
                               graph3JSON=graph3JSON, graph4JSON=graph4JSON)


    ############################################################################################################# Tajimas D and Nucleotide Diversity graphs########################################################
    if len(TD) != 0 and len(ND) != 0:

        # creating Tajimas D graph individual
        fig_1 = go.Figure()
        buttons = []
        i = 0

        # iterating through dict
        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_1.add_trace(obj)

            # args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(TD)
            args[i] = True

            # creating button object for each population
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

            # add button to list
            buttons.append(button)
            i += 1

        fig_1.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])
        # adding axis names
        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

        # graph for all populations Tajima D
        fig_2 = go.Figure()

        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)

        fig_2.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

        # Creating nucleotide diversity graph
        fig_3 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_3.add_trace(obj)
            args = [False] * len(ND)
            args[i] = True
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])
            buttons.append(button)
            i += 1

        fig_3.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_3.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Nucleotide Diversity")

        # graph for overall
        fig_4 = go.Figure()
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_4.add_trace(obj)
            fig_4.update_layout(
                xaxis_title='Position (Base pairs)',
                yaxis_title="Nucleotide Diversity")

        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)
        graph3JSON = json.dumps(fig_3, cls=plotly.utils.PlotlyJSONEncoder)
        graph4JSON = json.dumps(fig_4, cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('stats_2.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON,
                               graph3JSON=graph3JSON, graph4JSON=graph4JSON)


    ############################################################################################################ Haploid Diversity and Nucleotide Diversity graphs#########################################
    if len(HD) != 0 and len(ND) != 0:
        # Creating Halploid Diversity graph individual

        fig_1 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary
        for x in HD.items():
            # adding scatter line for each pop
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_1.add_trace(obj)

            # args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(HD)
            args[i] = True

            # create button object for each pop
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

            # add button to list
            buttons.append(button)
            i += 1

        fig_1.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

        # graph for overall Haploid Diversity_
        fig_2 = go.Figure()
        for x in HD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)
        fig_2.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

        # Creating nucleotide Diversity graph individual

        fig_3 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_3.add_trace(obj)
            args = [False] * len(ND)
            args[i] = True
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])
            buttons.append(button)
            i += 1

        fig_3.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_3.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Nucleotide Diversity")

        # graph for overall _

        fig_4 = go.Figure()
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_4.add_trace(obj)
            fig_4.update_layout(
                xaxis_title='Position (Base pairs)',
                yaxis_title="Nucleotide Diversity")



    # coverting to json objects
        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)
        graph3JSON = json.dumps(fig_3, cls=plotly.utils.PlotlyJSONEncoder)
        graph4JSON = json.dumps(fig_4, cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('stats_2.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON,
                                   graph3JSON=graph3JSON, graph4JSON=graph4JSON)




############################################################################################################## individual stats graphs #################################################################################################################

    # Only Tajimas D graph
    if len(TD) != 0:

        fig_1 = go.Figure()
        buttons = []
        i = 0

        #iterating through dict
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

        fig_1.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])
        # adding axis names
        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

    # graph for all populations Tajima D
        fig_2 = go.Figure()

        for x in TD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)

        fig_2.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Tajima's D")

        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('stats_3.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON)

    # Only Haploid Diversity graph
    if len(HD) != 0:
        fig_1 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary
        for x in HD.items():
            # adding scatter line for each pop
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_1.add_trace(obj)

            # args is a list of booleans that tells the buttons which trace to show on click
            args = [False] * len(HD)
            args[i] = True

            # create button object for each pop
            button = dict(label=x[0],
                          method="update",
                          args=[{"visible": args}])

            # add button to list
            buttons.append(button)
            i += 1

        fig_1.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

        # graph for overall Haploid Diversity_
        fig_2 = go.Figure()
        for x in HD.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)
        fig_2.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Haploid Diversity")

        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('stats_3.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON)

    # only Nucleotide Diversity graph

    else:
        fig_1 = go.Figure()
        buttons = []
        i = 0

        # iterating through dictionary
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_1.add_trace(obj)
            args = [False] * len(ND)
            args[i] = True

            button = dict(label=x[0],
                  method="update",
                  args=[{"visible": args}])

            buttons.append(button)
            i += 1

        fig_1.update_layout(
            updatemenus=[
                dict(
                    type="dropdown",
                    direction="down",
                    buttons=buttons)
            ])

        fig_1.update_layout(
            xaxis_title='Position (Base pairs)',
            yaxis_title="Nucleotide Diversity")

        # graph for overall _
        fig_2 = go.Figure()
        for x in ND.items():
            obj = go.Scatter(name=x[0], x=position, y=x[1])
            fig_2.add_trace(obj)
            fig_2.update_layout(
                xaxis_title='Position (Base pairs)',
                yaxis_title="Nucleotide Diversity")

        graph1JSON = json.dumps(fig_1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig_2, cls=plotly.utils.PlotlyJSONEncoder)

        return render_template('stats_3.html', graph1JSON=graph1JSON, graph2JSON=graph2JSON)

if __name__ == '__main__':
	app.run(debug=True)