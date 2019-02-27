
cat <<__SCRIPT__ > cluster.py
#!/bin/env python

import pandas as pd
import numpy as np
import seaborn

import fastcluster

z = pd.read_csv("$1", delimiter = "\t", header = None, skiprows = 1)

u = z.pivot(index=0, columns=1, values=2)

for i in np.arange(u.shape[0]):
    for j in np.arange(u.shape[1]):
        if(pd.isnull(u.iloc[i, j])):
            u.iloc[i, j] = u.iloc[j, i]

u[np.isnan(u)] = 0

g = seaborn.clustermap(u, row_cluster=True, col_cluster=True, xticklabels=False, yticklabels=False, cmap='Reds', vmin=0, vmax=1)
g.savefig("$2.png")

c = g.dendrogram_col.reordered_ind
r = g.dendrogram_row.reordered_ind

c = u.columns[c]
r = u.index[r]

with open("$2rows.txt", 'w') as fi:
    for item in r:
        fi.write("%s\n" % item)

with open("$2cols.txt", 'w') as fi:
    for item in c:
        fi.write("%s\n" % item)

__SCRIPT__
