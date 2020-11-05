import pandas as pd
import numpy as np
import math
import csv

from simulator import *
from simulator_plotting import *

pairs = [('CTX', 'SAM'), ('AM', 'AMC'), ('ZOX', 'CXM')]
pairs = [[dataset2.loc[item] for item in pair] for pair in pairs]

def results_to_csv(results, name):
    tree_file = 'march/{}_tree.csv'.format(name)
    pop_file = 'march/{}_pop.csv'.format(name)
    tr = results['trace']
    phylo = results['phylogeny']
    t_parents=[]
    t_idents=[]
    for ident, parent in enumerate(phylo):
        if parent != -2 and parent != -1:
            t_parents.append(parent)
            t_idents.append(ident)
    t_idents_ = ["{0:0>4b}".format(x) for x in t_idents]
    t_parents = ["{0:0>4b}".format(x) for x in t_parents]
    treedf = pd.DataFrame(data={'Parent': t_parents, 'Identity': t_idents_})
    treedf.to_csv(tree_file, quoting=csv.QUOTE_NONNUMERIC, index_label=False)

    sz = len(tr['0000'])
    df = pd.DataFrame(tr)
    gens=[]
    idents=[]
    pops=[]
    for i in range(sz):
        if i % 10 != 0: # only include every 10th timestep for smoother plot
            continue
        for ident, pop in enumerate(df.loc[i]):
            if ident in t_idents or ident == 0:
                gens.append(i)
                idents.append(ident)
                pops.append(math.log10(pop) if pop > 1 else 0) # log transform for more reasonable plot
    idents = ["{0:0>4b}".format(x) for x in idents]
    newdf = pd.DataFrame(data={'Generation': gens, 'Identity': idents, 'Population': pops})
    newdf.set_index('Generation', inplace=True)
    newdf.to_csv(pop_file, quoting=csv.QUOTE_NONNUMERIC)
    return tree_file, pop_file
results_to_csv(simulate(dataset1[2].iloc[4].tolist()), 'not_AMC')
