import sys, pprint
import numpy as np
import matplotlib.pyplot as plt

from io import StringIO   # StringIO behaves like a file object



dataDir = sys.argv[1]
#dataDir = './output1DFolder'
#pprint.pprint( dataDir)




#pprint.pprint( data.size )
#pprint.pprint( data.shape )
#pprint.pprint( data.dtype )
zmax = 50.0 #np.amax(data)
zmin = -100.0 #np.amin(data)
levels = np.linspace( zmin, zmax, num=151 )
cmap = plt.cm.jet

#
#figure;
# figure;
data = np.loadtxt( dataDir + '/ap10.dat');
fig1, axs1 = plt.subplots( constrained_layout=True )
CS1 = axs1.contourf( data , levels ,cmap = cmap )
axs1.set_axis_off()
cbar1 = fig1.colorbar( CS1 )
fig1.savefig( dataDir + '/apt10.png')

data = np.loadtxt( dataDir + '/ap20.dat');
fig2, axs2 = plt.subplots( constrained_layout=True )
CS2 = axs2.contourf( data , levels ,cmap = cmap )
axs2.set_axis_off()
cbar2 = fig2.colorbar( CS2 )
fig2.savefig( dataDir + '/apt20.png')




data = np.loadtxt( dataDir + '/ap30.dat');
fig3, axs3 = plt.subplots( constrained_layout=True )
CS3 = axs3.contourf( data , levels ,cmap = cmap )
axs3.set_axis_off()
cbar3 = fig3.colorbar( CS3 )
fig3.savefig( dataDir + '/apt30.png')




data = np.loadtxt( dataDir + '/ap100.dat');
fig4, axs4 = plt.subplots( constrained_layout=True )
CS4 = axs4.contourf( data , levels ,cmap = cmap )
axs4.set_axis_off()
cbar4 = fig4.colorbar( CS4 )
fig4.savefig( dataDir + '/apt100.png')




data = np.loadtxt( dataDir + '/ap200.dat');
fig5, axs5 = plt.subplots( constrained_layout=True )
CS5 = axs5.contourf( data , levels ,cmap = cmap )
axs5.set_axis_off()
cbar5 = fig5.colorbar( CS5 )
fig5.savefig( dataDir + '/apt200.png')



data = np.loadtxt( dataDir + '/ap300.dat');
fig6, axs6 = plt.subplots( constrained_layout=True )
CS6 = axs6.contourf( data , levels ,cmap = cmap )
axs6.set_axis_off()
cbar6 = fig6.colorbar( CS6 )
fig6.savefig( dataDir + '/apt300.png')

#plt.show()
