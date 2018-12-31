import sys, pprint
import numpy as np
import matplotlib.pyplot as plt

from io import StringIO   # StringIO behaves like a file object

dataDir = sys.argv[1]
#dataDir = './output0DFolder'
#pprint.pprint( dataDir)

fname = dataDir + '/vm_1Hz.txt'
#pprint.pprint(fname)

ty = np.loadtxt(fname) #'vm_1Hz.txt');
#pprint.pprint(ty)

ts = ty[:,0]/1e3

tss = ts
#%%
y = ty[:,1:];

syids = 8;
yids = np.array( [29, 30, 31, 32, 33, 35, 36, 37, 38] )

ynid = np.zeros(217, dtype=int );


for id in range(0, (syids+1) ):
    ynid[ yids[id] ] = id

#end


mycolor = 'b'

ar = np.loadtxt(dataDir + '/allresult_1Hz.txt')

tArray = ar[:,0]
INa = ar[:,1]
I_Ca_store = ar[:,2]
Ikur1 = ar[:,3]
Ikur2 = ar[:,4]
Iss = ar[:,5]
Jserca = ar[:,6]
Jleak = ar[:,7:9]
Nai = ar[:,9]
Incx = ar[:,10]


#%%
#% Membrane potential
#figure;
fig1, axs1 = plt.subplots(figsize=(12.8,9.6))
axs1.plot( ts, y[:, ynid[38] ], mycolor)
#hold on;
axs1.set_xlabel('Time (sec)', fontsize=20)
axs1.set_ylabel('Membrane Potential (mV)', fontsize=20);

#%saveas(gcf, 'vm', 'fig');
#saveas(gcf, [fullPath, 'vm'], 'png');
# fig1.set_figheight(14.4)
# fig1.set_figwidth(19.2)
fig1.savefig( dataDir + '/vm.png' )

#%%
#% LCC current (ICa)
#figure;
fig2, axs2 = plt.subplots(figsize=(12.8,9.6))
axs2.plot( tss, I_Ca_store, mycolor)
axs2.set_xlabel('Time (sec)', fontsize=20)
axs2.set_ylabel('I$_{Ca}$ (pA/pF)', fontsize=20)

#%saveas(gcf, 'ICa', 'fig');
#saveas(gcf, [fullPath, 'ICa'], 'png');
fig2.savefig( dataDir + '/ICa.png' )



#%%
#% CaSRT & Caj
#figure;
fig3, axs3 = plt.subplots(1, 3, figsize=(12.8, 9.6) )
fig3.subplots_adjust(wspace=0.5)
#subplot(1,3,1)
axs3[0].plot( ts, ( y[:,ynid[29]] + y[:, ynid[30]] ), mycolor )
axs3[0].set_xlabel('Time (sec)', fontsize=20)
axs3[0].set_ylabel('[Ca]$_{SRT}$ (mM)', fontsize=20)
#subplot(1,3,2)
axs3[1].plot( ts, y[:, ynid[35] ]*1e3, mycolor)
axs3[1].set_xlabel('Time (sec)', fontsize=20)
axs3[1].set_ylabel('Ca Dyad ($\mu$M)', fontsize=20)
#subplot(1,3,3)
axs3[2].plot( ts, y[:, ynid[36] ], mycolor)
axs3[2].set_xlabel('Time (sec)', fontsize=20)
axs3[2].set_ylabel('Ca sl (mM)', fontsize=20)
#
#%saveas(gcf, 'CaSRT_Caj', 'fig');
#saveas(gcf, [fullPath, 'CaSRT_Caj'], 'png');
fig3.savefig( dataDir + '/CaSRT_Caj.png' )

#
#%%
#% Cai
#figure;
fig4, axs4 = plt.subplots( figsize=(12.8, 9.6) )
axs4.plot( ts, y[:, ynid[37] ], mycolor )
axs4.set_xlabel('Time (sec)', fontsize=20 )
axs4.set_ylabel('[Ca]$_i$ ($\mu$M)', fontsize=20 )
#
#%saveas(gcf, 'Cai', 'fig');
#saveas(gcf, [fullPath, 'Cai'], 'png');
fig4.savefig( dataDir + '/Cai.png' )



#% Ikur
#figure;
fig5, axs5 = plt.subplots( figsize=(12.8, 9.6) )
axs5.plot( tss, Ikur1, mycolor )
axs5.set_xlabel('Time (sec)', fontsize=20 )
axs5.set_ylabel('I$_{ku1}$ (pA/pF)', fontsize=20 )
#
fig5.savefig( dataDir + '/Iku1.png' )
# saveas(gcf, [fullPath, 'Iku1'], 'png');

#% Ikur
#figure;
fig6, axs6 = plt.subplots( figsize=(12.8, 9.6) )
axs6.plot( tss, Ikur2, mycolor )
axs6.set_xlabel('Time (sec)', fontsize=20 )
axs6.set_ylabel('I$_{ku2}$ (pA/pF)', fontsize=20 )
#
fig6.savefig( dataDir + '/Iku2.png' )

#
#% INa
#figure;
fig7, axs7 = plt.subplots( figsize=(12.8, 9.6) )
axs7.plot( tss, INa, mycolor )
axs7.set_xlabel('Time (sec)', fontsize=20 )
axs7.set_ylabel('I$_{Na}$ (pA/pF)', fontsize=20 );
#
#%saveas(gcf, 'INa', 'fig');
#saveas(gcf, [fullPath, 'INa'], 'png');
fig7.savefig( dataDir + '/INa.png' )



#% Iss
#figure;
fig8, axs8 = plt.subplots( figsize=(12.8, 9.6) )
axs8.plot( tss, Iss, mycolor )
axs8.set_xlabel('Time (sec)', fontsize=20 )
axs8.set_ylabel('I$_{ss}$ (pA/pF)', fontsize=20 );
#
fig8.savefig( dataDir + '/Iss.png' )




#% [Na]
#figure;
fig9, axs9 = plt.subplots( figsize=(12.8, 9.6) )
axs9.plot( tss, Nai, mycolor )
axs9.set_xlabel('Time (sec)', fontsize=20 )
axs9.set_ylabel('[Na] (mmol/L relevant compartment', fontsize=20 );
#saveas(gcf, [fullPath, 'Na'], 'png');
fig9.savefig( dataDir + '/Na.png' )


#% I_NCX
#figure;
fig10, axs10 = plt.subplots( figsize=(12.8, 9.6) )
axs10.plot( tss, Incx, mycolor) # %legend('I_N_C_X');
axs10.set_xlabel('Time (sec)', fontsize=20 )
axs10.set_ylabel('I$_{NCX}$ (pA/pF)', fontsize=20 )

#%saveas(gcf, 'I_NCX', 'fig');
#saveas(gcf, [fullPath, 'I_NCX'], 'png');
fig10.savefig( dataDir + '/I_NCX.png' )




