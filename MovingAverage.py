import math
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import glob
import os

coordinate_filepath = '/Users/shohtatakami/physics/COMETDS/LaserTracker/Interpolated_Map/X_0Y_0.csv_Map.csv'
column_names = ['encoder', 'x', 'y', 'z', 'Bx', 'By', 'Bz', 'B']
coordinate = pd.read_csv(coordinate_filepath, sep=',',header=None, names = column_names)
x = coordinate['x']
y = coordinate['y']
z = coordinate['z']
xzangle = []
yzangle = []
for i in range(1,len(z),1):
    rail_thetaXZ = math.atan((x[i]-x[i-1])/(z[i]-z[i-1]))*(180/math.pi)
    rail_thetaYZ = math.atan((y[i]-y[i-1])/(z[i]-z[i-1]))*(180/math.pi)
    print(rail_thetaXZ)
    xzangle.append(rail_thetaXZ)
    yzangle.append(rail_thetaYZ)

z_for_plot = z[:-1]

""""
for i in range(0,len(z),1):
    print(f"x : {x[i]} y : {y[i]} z : {z[i]}")
"""

    


#---------------------------------------------------------------------
#移動平均
#---------------------------------------------------------------------

def moving_average(value,term):
    """
    data:平滑化するデータ
    term:移動平均を取る区間
    """
    kernel = np.ones(term)/float(term)
    return np.convolve(value, kernel, mode='valid')

term = 10
xzangle_ma = moving_average(xzangle, term)
yzangle_ma = moving_average(yzangle, term)
z_ma = moving_average(z, term+1)#term+1?

print(len(x))
print(len(y))
print(len(z))
print(len(z_for_plot))
print(len(z_ma))
print(len(xzangle))
print(len(yzangle))
print(len(xzangle_ma))
print(len(yzangle_ma))

fig, ax = plt.subplots(2,1,figsize=(12,10))
#z vs x
ax[0].scatter(z_for_plot,xzangle,marker='o',color='red',label='original')
ax[0].scatter(z_ma,xzangle_ma,color='blue',label='moving avraged')
ax[0].set_title('XZ Angle')
ax[0].set_xlabel('z (mm)')
ax[0].set_ylabel('xz angle (°)')
ax[0].legend()
ax[0].grid(True)

#z vs y
ax[1].scatter(z_for_plot,yzangle,marker='o',color='red',label='original')
ax[1].scatter(z_ma,yzangle_ma,color='blue',label='moving avraged')
ax[1].set_title('YZ Angle')
ax[1].set_xlabel('z (mm)')
ax[1].set_ylabel('yz angle (°)')
ax[1].legend()
ax[1].grid(True)

plt.tight_layout
plt.show()
    