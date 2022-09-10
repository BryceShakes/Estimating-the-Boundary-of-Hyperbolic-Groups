# %%

# libraries
import numpy as np
import pandas as pd
from sklearn import manifold
import matplotlib.pyplot as plt
import os

# %%        my function that should work


# takes the input of an n x n array (symmetrical) with pairwise distances, called a distance or dissimilarity matrix
# performs MDS, either metric or non based on input, with fixed (although tuned) parameters
# returns the points as cartesian co-ords of specified dimension in Euclidean space
# also plots the graph if dimensions = 2
def plot_from_matrix(dist, metric=True, title = None, dimensions = 2):
    if metric == False:
        mds = manifold.MDS(n_components=dimensions,
                       max_iter=100000,
                       n_init = 5,
                       eps=1e-15,
                       metric = False,
                       dissimilarity="precomputed")
        new = mds.fit(dist).embedding_
    elif metric:
        mds = manifold.MDS(n_components=dimensions,
               max_iter=10000,
               n_init = 5,
               eps=1e-3,
               metric = True,
               dissimilarity="precomputed",
               verbose = 1)
        new = mds.fit(dist).embedding_
    if dimensions == 2:
        plt.scatter(new[:,0], new[:,1])
        if title != None:
            plt.title(title)
        plt.show()
    return new

# takes the CSV output from distance-matrix-generation-functions.g and converts into a numpy array
def csv_to_matrix(path):
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

# combines plot_from_matrix and csv_to_matrix. 
# supply with path as location of CSV file from GAP and equal parameters to 
def plot_from_csv(path, metric = True, title = None,  dimensions =2 ):
    arr = csv_to_matrix(path)
    return plot_from_matrix(arr, metric=metric, title=title, dimensions = dimensions)


# not my code, edited slightly to nest in function from https://stackoverflow.com/questions/54055761/checking-triangle-inequality-in-a-massive-numpy-matrix
def check_triangle_ineq(D):
    for i in range(len(D)):
        for j in range(i):
            if not all(D[i,j] <= D[i,:] + D[:,j]):
                return False
    return True

# count the number of triangles in of D that violate the traingle inequality. from a total of nC3, where D is an nxn matrix
# again not my code, was altered from the check_triangle_ineq
def measure_triangle_ineq(D):
    N = len(D)
    test = 0
    for i in range(N):
        for j in range(i):
            test += sum(D[i,j] > D[i,:] + D[:,j])
    return test

# %% genus 2 surface

path = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/genus_2_vis.csv"

plot_from_csv(path, False);
    

# %%

#  C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1
path = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/Surface_500.csv"
path2 = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/free_vis_mat_noid.csv" 


a1 = csv_to_matrix(path)
a2 = csv_to_matrix(path2)

plot_from_matrix(a1, metric = True)
plot_from_matrix(a2, metric = True, title = "Free group, n = 1000, epsilon = 1")

#bryce_draw_from_raw(path)


# %% testing speeds and tuning manifold params

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
    print(check_triangle_ineq(csv_to_matrix(path)));

for i in list:
    path = r"C:/Users/Bryce/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc/LocalState/rootfs/home/bryce/gap-4.11.1/"+i+".csv"
    plot_from_csv(path, True, title = i);

D = np.array([[0, 3, 5, 9, 2],
             [1, 0, 4, 7, 1 ],
             [5, 4, 0, 8, 6],
             [9, 7, 8, 0, 9],
             [2, 1, 6, 9, 0]])

P = np.array([[0, 4, 100],
             [4, 0, 5,],
             [100, 5, 0,]])

measure_triangle_ineq(a1)

#%%

import decimal

a = decimal.Decimal(4*3**999)
format(a, '.2e')

b = decimal.Decimal(1e+120)

c = decimal.Decimal((1000*8*4*3**999) /b )
format(c, '.2e')
