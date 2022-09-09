# %%

from time import time
import numpy as np
import pandas as pd
from sklearn import manifold
import matplotlib.pyplot as plt
import os

# %%        my function that should work


# takes the input of an n x n array (symmetrical) with pairwise distances

def Bryce_picture(dist, metric=True, title = None):
    if metric == False:
        mds = manifold.MDS(n_components=2,
                       max_iter=100000,
                       n_init = 5,
                       eps=1e-15,
                       metric = False,
                       dissimilarity="precomputed")
        new = mds.fit(dist).embedding_
    elif metric:
        mds = manifold.MDS(n_components=2,
               max_iter=10000,
               n_init = 5,
               eps=1e-3,
               metric = True,
               dissimilarity="precomputed",
               verbose = 1)
        new = mds.fit(dist).embedding_
    plt.scatter(new[:,0], new[:,1])
    if title != None:
        plt.title(title)
    plt.show()

def csv_to_mat(path):
    assert os.path.isfile(path)
    df = pd.read_csv(path)
    if df.shape[1] > 1:
        arr = []
        for i in range(df.shape[1]):
            arr.append(eval(df.iloc[0,i]))
        arr = np.array(arr)
    elif df.shape[1] == 1:
        arr = np.array(eval(df.iloc[0,0]))
    return arr

def bryce_draw_from_raw(path, metric = True, title = None):
    arr = csv_to_mat(path)
    Bryce_picture(arr, metric=metric, title=title)

# not my code, altered ever so slightly from https://gist.github.com/tuelwer/b7ad6d2e69a823d3302f2f9996f783e6
def check_triangle_inequality_verbose(D):
    """ Returns true iff the matrix D fulfills
    the triangle inequaltiy.
    """
    n = len(D)
    valid = True
    for i in range(n):
        for j in range(i, n):
            for k in range(n):
                if k == j or k == i:
                    continue
                if D[i][j] > D[i][k] + D[k][j]:
                    print('Invalid triple:', D[i][j],  D[i][k], D[k][j])
                    valid = False
    return valid

# not my code, edited slightly to nest in function from https://stackoverflow.com/questions/54055761/checking-triangle-inequality-in-a-massive-numpy-matrix
def check_triangle_inequality(D):
    N = len(D)
    test = True
    for i in range(N):
        for j in range(i):
            test = test & all(D[i,j] <= D[i,:] + D[:,j])
            if not test :
                return test
    return test

# count the number of triangles in of D that violate the traingle inequality. from a total of nC3, where D is an nxn matrix
def measure_triangle_ineq(D):
    N = len(D)
    test = 0
    for i in range(N):
        for j in range(i):
            test += sum(D[i,j] > D[i,:] + D[:,j])
    return test

# %% genus 2 surface

path = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/genus_2_vis.csv"

bryce_draw_from_raw(path, True);
    

# %%

#  C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1
path = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/Surface_500.csv"
path2 = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/free_vis_mat_noid.csv" 


a1 = csv_to_mat(path)
a2 = csv_to_mat(path2)

Bryce_picture(a1, metric = True)
Bryce_picture(a2, metric = True, title = "Free group, n = 1000, epsilon = 1")

#bryce_draw_from_raw(path)


# %% testing speeds

testing = [
    [False, 1e-15, 100000],
    [False, 1e-9, 100000],
    [False, 1e-15, 1000],
    [False, 1e-9, 1000],
    [True, 1e-15, 100000],
    [True, 1e-9, 100000],
    [True, 1e-3, 100000],
    [True, 1e-15, 10000],
    [True, 1e-9, 10000],
    [True, 1e-3, 10000],
    [True, 1e-15, 1000],
    [True, 1e-9, 1000],
    [True, 1e-3, 1000],
    ]

now = time()
time_arr = []
for i in testing :
    mds = manifold.MDS(n_components=2,
           max_iter = i[2],
           n_init = 5,
           eps = i[1],
           metric = i[0],
           dissimilarity = "precomputed")
    new = mds.fit(a1).embedding_
    plt.scatter(new[:,0], new[:,1])
    plt.title(str(i))
    plt.show()
    time_arr.append([i,time() - now - sum([x[1] for x in time_arr])])
    
# %% testing with genus 2 surface groups

list = ["surface_l_100_e_05", "surface_l_100_e_1", "surface_l_100_e_2","surface_l_300_e_05","surface_l_300_e_1",
       "surface_l_300_e_2","surface_l_500_e_05","surface_l_500_e_1","surface_l_500_e_2"]

for i in list:
    path = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/"+i+".csv"
    print(check_triangle_inequality(csv_to_mat(path)));

for i in list:
    path = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/"+i+".csv"
    bryce_draw_from_raw(path, True, title = i);

D = np.array([[0, 3, 5, 9, 2],
             [1, 0, 4, 7, 1 ],
             [5, 4, 0, 8, 6],
             [9, 7, 8, 0, 9],
             [2, 1, 6, 9, 0]])

P = np.array([[0, 4, 100],
             [4, 0, 5,],
             [100, 5, 0,]])


    
measure_triangle_ineq(a1)
