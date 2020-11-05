import pandas as pd
import numpy as np
#import igraph as ig

from plotly.offline import plot, iplot
import plotly.graph_objs as go

def get_dataset1():
    dataset = []
    for i in [1,2,3]:
        df = pd.read_csv(f'../data/table{i}.csv', dtype={'Genotype': str}).set_index('Genotype').T
        df = df.rename_axis('Label')
        dataset.append(df)
    return dataset

dataset1 = get_dataset1()

def get_dataset2():
    df = pd.read_csv('../data/barlow.csv', index_col=0)
    df = df.sort_index(axis=1)
    df = df.rename_axis('Genotype', axis='columns')
    return df

dataset2 = get_dataset2()

def genotype_strings(num):
    return(["{0:0>4b}".format(n) for n in range(num)])

def plot_simulation(results, genotypes_to_plot=[], carrying_capacity=9):
    trace = results['trace']
    if genotypes_to_plot and all(isinstance(x, int) for x in genotypes_to_plot):
        genotypes_to_plot = ["{0:0>4b}".format(n) for n in genotypes_to_plot]
    if not genotypes_to_plot:
        genotypes_to_plot = list(
            set(results['greedy_path']) | set(results['actual_path'])
        )
        genotypes_to_plot.sort(key=lambda g: np.argmax(trace[g] > 0))
    y_max = carrying_capacity + 1
    timesteps = len(next(iter(trace.values())))
    #           length of an arbitrary entry from dictionary
    index = np.array(range(0, timesteps), np.int64)
    
    # data
    data = []
    for g in genotypes_to_plot:
        data.append(go.Scatter(
            x=index,
            y=trace[g],
            name=g,
            line=dict(width=1)
        ))

    # vertical lines
    for crit in ['T_1', 'T_d', 'T_f']:
        if results[crit] != -1:
            data.append(go.Scatter(
                x=[results[crit], results[crit]],
                y=[0, 10**y_max],
                line=dict(color='black', width=1),
                name = crit,
                mode = 'lines',
                showlegend = False
            ))

    # vertical line annotations
    vlines = []
    for crit, text in zip(['T_1', 'T_d', 'T_f'],
                          ['<i>T</i><sub>1</sub>',
                           '<i>T</i><sub>d</sub>',
                           '<i>T</i><sub>f</sub>']):
        if results[crit] != -1:
            vlines.append(dict(
                x=results[crit],
                y=y_max*.975,
                xref='x',
                yref='y',
                text=text,
                showarrow=False,
                xanchor='left'
            ))
        
    # layout
    layout = go.Layout(
        xaxis=dict(
            title='Timestep'
        ),
        yaxis=dict(
            type='log',
            range=[0,y_max],
            exponentformat='power',
            title='Abundance of each genotype',
            tickfont=dict(size=10)
        ),
        width=750, #600
        height=500, #400
        margin=dict(t=30, b=50),
        annotations=vlines
    )
        
    fig = go.Figure(data=data, layout=layout)
    
    iplot(fig, show_link=False)

# an example:
#landscape = dataset1[2].loc['53.60μM'].tolist()
#results = simulate(landscape, timesteps=2000)
#plot_simulation(results, ["0000", "0010", "0110", "1110"])

def plot_path(path, num_genotypes):
    g = ig.Graph()
    g.add_vertices(num_genotypes)
    g.vs["name"] = genotype_strings(num_genotypes)

    for gen in range(num_genotypes):
        neighbors = [gen ^ 1 << i for i in range(4)]
        for n in neighbors:
            if not g.are_connected(gen, n):
                g.add_edge(gen, n)

    g.es["on_path"] = False
    for i, gen in enumerate(path[:-1]):
        e = g.get_eid(g.vs.find(gen), g.vs.find(path[i+1]))
        g.es[e]["on_path"] = True

    visual_style = {
        'layout': g.layout("rt"),
        'vertex_label': g.vs["name"],
        'bbox': (300,300),
        'vertex_size': 30,
        'vertex_label_size': 10,
        'vertex_color': "Powder Blue",
        'edge_width': [1 + 2 * int(on_path) for on_path in g.es["on_path"]]
    }
    layout = g.layout("rt")
    return ig.plot(g, **visual_style)
