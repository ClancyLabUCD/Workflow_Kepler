import sys, pprint
import numpy as np
import matplotlib.pyplot as plt

from io import StringIO   # StringIO behaves like a file object

dataDir = sys.argv[1]
#dataDir = './output1DFolder'
#pprint.pprint( dataDir)



# figure;
a = np.loadtxt( dataDir + '/ECGs.txt');

fig1, axs1 = plt.subplots(figsize=(12.8,9.6))

axs1.plot(a[:,0],a[:,1],'k');

#saveas(gcf, 'ECG', 'fig');
#saveas(gcf, [fullPath, 'ECG'], 'png');
fig1.savefig( dataDir + '/ECG.png')
