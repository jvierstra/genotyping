

cat <<__SCRIPT__ > cluster.py
#!/bin/env python

import pandas as pd
import numpy as np

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, fcluster
import fastcluster

z = pd.read_csv("results/relatedness.dsnums.txt", delimiter = "\t", header = None)

u = z.pivot(index=0, columns=1, values=2)

for i in np.arange(u.shape[0]):
    for j in np.arange(u.shape[1]):
        if(np.isnan(u.iloc[i, j])):
            u.iloc[i, j] = u.iloc[j, i]

u[np.isnan(u)] = 0

D = scipy.spatial.distance.pdist(u, 'correlation')
L = fastcluster.linkage(D, method = 'complete')
clusters = fcluster(L, .5, 'distance')

#Write file

#Diagnostic plots


dendro = dendrogram(L, no_plot = False)
o2 = dendro["leaves"]

clusters = fcluster(L, .5, 'distance')
__SCRIPT__