# %%

import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd

NUMBER_OF_KMERS = 136
NUMBER_OF_SEQ = 238
NUMBER_OF_CLUSTERS = 10
SEQ_IN_CLUSTER = [23,25,24,22,24,24,24,24,25,23]

X = np.empty((0,NUMBER_OF_KMERS), int)
temp = []
#print(X)
for i in range(1,NUMBER_OF_SEQ+1):
    filename = 'mer_counts_dumps_' + str(i) +'.fa'
    with open(filename,'r') as file:
        for number, line in enumerate(file):
            if number % 2 == 0:
                temp.append(int(line[1:]))
    X = np.append(X, np.array([temp]), axis=0)
    temp.clear()

db = DBSCAN(eps = 34).fit(X) #empirical value for eps

print(db.labels_)

pca = PCA(n_components=2)
pca.fit(X)
a = pca.fit_transform(X)

df = pd.DataFrame({'x': a[:,0],
                   'y': a[:,1],
                   'z': db.labels_})

print(type(db.labels_))
truelabels = []

for template in range(1,11):
    range_up = sum(SEQ_IN_CLUSTER[0:(template-1)]) + 1
    range_down = sum(SEQ_IN_CLUSTER[0:template]) + 1
    for i in range(range_up,range_down):
        truelabels.append(template)

groups = df.groupby('z')
#color codes for 10 clusters
c = ['k','orange','lime','r','darkorchid','peru','magenta','grey','greenyellow','aqua','royalblue']

color_ind = 0
for name, group in groups:
    plt.plot(group.x, group.y, marker='o', linestyle='', markersize=5, label=name, color=c[color_ind])
    color_ind = color_ind + 1

plt.legend()
plt.show()

df1 = pd.DataFrame({'x': a[:,0],
                   'y': a[:,1],
                   'z': truelabels})
print(df1)
groups = df1.groupby('z')
#color codes for 10 clusters
c = ['orange','lime','r','darkorchid','peru','magenta','grey','greenyellow','aqua','royalblue']

color_ind = 0
print(len(groups))
for name, group in groups:
    plt.plot(group.x, group.y, marker='o', linestyle='', markersize=5, label=name, color=c[color_ind])
    color_ind = color_ind + 1

plt.legend()
plt.show()
# %%
