from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics import silhouette_score
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import cosine_similarity

truelist = []
predlist = []

def purity():
    score = 0
    total = 0
    with open('DBoutput.txt','r') as file1, open('trueclusters.txt','r') as file2:
        for line1, line2 in zip(file1, file2):
            total = total + 1
            if(line1==line2):
                score = score + 1
    return score/total

def NMI():
    score  = normalized_mutual_info_score(truelist, predlist)
    return score

def silhouette(X, intra, inter):
    score  = silhouette_score(X, predlist)
    alt_score = (intra-inter)/max(intra,inter)
    print('alternate silhoutte score: ', alt_score)
    return score

def intracluster_similarity(all_list):

    cluster_indexes = []
    cluster_items = []
    all_sum = 0
    for i in range(0,10):
        # print('when i is ',i)
        cluster_indexes = [ind for ind, label in enumerate(predlist) if label == i]
        # print(cluster_indexes)
        cluster_items = [all_list[ind] for ind in cluster_indexes]
        # print('cluster items len ',len(cluster_items))
        # print('cluster items len ',len(cluster_items[0]))
        ed = cosine_similarity(cluster_items,cluster_items)
        
        # print(len(ed))
        # print(len(ed[0]))
        total = 0
        for j in ed:
            total = total + sum(j)

        avg_dist = total/(len(ed)*len(ed[0]))
        # print(total)
        # print(avg_dist)
        all_sum = all_sum + avg_dist
    return all_sum/10

def intercluster_similarity(all_list):
    all_sum = 0
    all_sum = 0
    cnt = 0
    for i in range(0,10):
        for j in range(0,10):
            if i==j:
                continue

            # print('when i is ',i,' and j is ',j)
            cluster_indexes1 = [ind for ind, label in enumerate(predlist) if label == i]
            cluster_indexes2 = [ind for ind, label in enumerate(predlist) if label == j]
            # print(cluster_indexes1)
            # print(cluster_indexes2)
            cluster_items1 = [all_list[ind] for ind in cluster_indexes1]
            cluster_items2 = [all_list[ind] for ind in cluster_indexes2]
            # print('cluster items 1 len ',len(cluster_items1))
            # print('cluster items 1 len ',len(cluster_items1[0]))
            ed = cosine_similarity(cluster_items1,cluster_items2)
            # print(len(ed))
            # print(len(ed[0]))
            total = 0
            for k in ed:
                total = total + sum(k)
            avg_dist = total/(len(ed)*len(ed[0]))
            # print(total)
            # print(avg_dist)
            all_sum = all_sum + avg_dist
            cnt = cnt + 1
    #print('count is ',cnt)
    return all_sum/cnt

with open('DBoutput.txt','r') as file1, open('trueclusters.txt','r') as file2:
    for line1, line2 in zip(file1, file2):
        truelist.append(int(line2))
        predlist.append(int(line1))

X = np.empty((0,136), int)
temp = [] #1d list
all_list = [] #2d list
for i in range(0,236):
    filename = 'mer_counts_dumps_' + str(i) +'.fa'
    with open(filename,'r') as file:
        for number, line in enumerate(file):
            if number % 2 == 0:
                temp.append(int(line[1:]))
        all_list.append(temp.copy())

    X = np.append(X, np.array([temp]), axis=0)
    temp.clear()

print("purity: ",purity())
print("NMI: ",NMI())
intra = intracluster_similarity(all_list)
inter = intercluster_similarity(all_list)
print("intra cluster similarity: ", intra)
print("inter cluster similarity: ", inter)
print("silhouette score: ",silhouette(X, intra, inter))