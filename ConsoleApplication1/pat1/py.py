import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
import os
import meshio
import csv
from pygem import RBFFactory
import copy
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl
def loadCV(filename):
    # initializing the titles and rows list
    fields = []
    rows = []

    # reading csv file
    with open(filename, 'r') as csvfile:
        # creating a csv reader object
        csvreader = csv.reader(csvfile)
        # extracting field names through first row
        fields = next(csvreader)
        # extracting each data row one by one
        for row in csvreader:
            rows.append([float(i) for i in row])
    return np.array(rows)
#colors = ["darkorange", "gold", "lawngreen", "lightseagreen"]
#cmap1 = LinearSegmentedColormap.from_list("mycmap", colors)
d = loadCV("C:/Users/polit/source/repos/ConsoleApplication1/ConsoleApplication1/pat1/distancePoints_1.csv")
dn = [];
z=0;
for i in d:
    if ((i[0]>-0.8) & (i[0]!=0)):
        dn.append(i[0])
        #print(i[0])
    else:
        z=z+1;
print(z)

# Displaying the graph
cmapM = plt.cm.winter
cmapMax = plt.cm.autumn

fig, axis = plt.subplots(figsize =(10, 5))
cnts, values, bars =axis.hist(dn, bins=[-0.173*3,-0.173*2,-0.173,0,0.173,0.173*2,0.173*3])
colors = ["blue", "aqua", "#47d147"]
cmap1 = LinearSegmentedColormap.from_list("mycmap", colors)
colors = ["#47d147", "yellow", "red"]
cmap = LinearSegmentedColormap.from_list("mycmap", colors)
for i, (cnt, value, bar) in enumerate(zip(cnts, values, bars)):
    #print(values)
    if value > 0:
        bar.set_facecolor(cmap(value/values.max()))
    if value < 0:
        bar.set_facecolor(cmap1(1-value/values.min()))
plt.title("Patient 1 Aneurysm Distnace between points for time step 6 and 10")
plt.xlabel("Movement (mm)")
plt.ylabel("Count")
plt.show()
