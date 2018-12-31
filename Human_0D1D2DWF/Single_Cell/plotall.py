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

ynid = np.zeros(206, dtype=int )

for id in range(0, (syids+1) ):
    ynid[ yids[id] ] = id

#end


mycolor = 'b';

ar = np.loadtxt(dataDir + '/allresult_1Hz.txt')

tArray = ar[:,0];
I_Ca_store = ar[:,1];
Ito = ar[:,2];
Itof = ar[:,3];
Itos = ar[:,4];
INa = ar[:,5];
IK1 = ar[:,6];
s1 = ar[:,7];
k1 = ar[:,8];
Jserca = ar[:,9];
Iks = ar[:,10];
Ikr = ar[:,11];
Jleak = ar[:,12:14];
ICFTR = ar[:,14];
Incx = ar[:,15];

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

#
#% Ito
#figure;
fig5, axs5 = plt.subplots( figsize=(12.8, 9.6) )
axs5.plot( tss, Ito , mycolor )
axs5.set_xlabel( 'Time (sec)', fontsize=20 )
axs5.set_ylabel( 'I$_{to}$ (pA/pF)', fontsize=20 )
#% legend('I_t_o','I_to_f','I_to_s');
#
#%saveas(gcf, 'Ito', 'fig');
#saveas(gcf, [fullPath, 'Ito'], 'png');
fig5.savefig( dataDir + '/Ito.png' )

#
#% INa 
#figure;
fig6, axs6 = plt.subplots( figsize=(12.8, 9.6) )
axs6.plot( tss, INa, mycolor )
axs6.set_xlabel('Time (sec)', fontsize=20 )
axs6.set_ylabel('I$_{Na}$ (pA/pF)', fontsize=20 );
#
#%saveas(gcf, 'INa', 'fig');
#saveas(gcf, [fullPath, 'INa'], 'png');
fig6.savefig( dataDir + '/INa.png' )

#
#
#%%
#%%
#
#% IKs and ICFTR
#figure;
fig7, axs7 = plt.subplots( 2, 1, figsize=(12.8, 9.6) )
fig7.subplots_adjust(hspace=0.5)
#subplot(2,1,1)
axs7[0].plot( tss, Iks )
# %legend('I_K_s (pA/pF)')
axs7[0].set_xlabel('Time (sec)', fontsize=20 )
axs7[0].set_ylabel('I$_{Ks}$ (pA/pF)', fontsize=20 );
#
#subplot(2,1,2);
axs7[1].plot( tss, ICFTR )
axs7[1].legend(['I$_{CFTR}$'] , fontsize=20)
# axs7[1].set_xlabel('Time (sec)', fontsize=20 )
# axs7[1].set_ylabel('I$_{Ks_{ICFTR}}$ (pA/pF)', fontsize=20 );

#%saveas(gcf, 'IKs_ICFTR', 'fig');
#saveas(gcf, [fullPath, 'IKs_ICFTR'], 'png');
fig7.savefig( dataDir + '/IKs_ICFTR.png' )

#
#
#% IKr and IK1
#figure;
fig8, axs8 = plt.subplots( 2, 1, figsize=(12.8, 9.6) )
fig8.subplots_adjust(hspace=0.5)
#subplot(2,1,1)
axs8[0].plot( tss, Ikr )
# %legend('I_K_r (pA/pF)')
axs8[0].set_xlabel('Time (sec)', fontsize=20 )
axs8[0].set_ylabel('I$_{Kr}$ (pA/pF)', fontsize=20 )
#
#subplot(2,1,2);
axs8[1].plot(tss,IK1)
# %legend('I_K_1 (pA/pF)')
axs8[1].set_xlabel('Time (sec)', fontsize=20 )
axs8[1].set_ylabel('I$_{K1}$ (pA/pF)', fontsize=20)
#
#
#%saveas(gcf, 'IKr_IK1', 'fig');
#saveas(gcf, [fullPath, 'IKr_IK1'], 'png');
fig8.savefig( dataDir + '/IKr_IK1.png' )

#
#
#% [Na]
#figure;
fig9, axs9 = plt.subplots( 1, 3, figsize=(12.8, 9.6) )
fig9.subplots_adjust(wspace=0.5)
#subplot(1,3,1);
axs9[0].plot(ts,y[:, ynid[31] ], mycolor )
axs9[0].legend(['[Na]$_j$'], fontsize=20 )
#subplot(1,3,2);
axs9[1].plot(ts,y[:, ynid[32] ], mycolor )
axs9[1].legend(['[Na]$_{sl}$'], fontsize=20 )
#subplot(1,3,3);
axs9[2].plot( ts, y[:, ynid[33] ], mycolor )
axs9[2].legend(['[Na]$_i$'], fontsize=20);
axs9[2].set_ylabel('[Na] (mmol/L) relevant compartment', fontsize=20 )
#
#%saveas(gcf, 'Na', 'fig');
#saveas(gcf, [fullPath, 'Na'], 'png');
fig9.savefig( dataDir + '/Na.png' )

#
#
#% I_NCX
#figure;
fig10, axs10 = plt.subplots( figsize=(12.8, 9.6) )
axs10.plot( tss, Incx, mycolor) # %legend('I_N_C_X');
axs10.set_xlabel('Time (sec)', fontsize=20 )
axs10.set_ylabel('I$_{NCX}$ (pA/pF)', fontsize=20 )

#%saveas(gcf, 'I_NCX', 'fig');
#saveas(gcf, [fullPath, 'I_NCX'], 'png');
fig10.savefig( dataDir + '/I_NCX.png' )

#
#
#% RyR fluxes
#figure;
fig11, axs11 = plt.subplots( 3, 1, figsize=(12.8, 9.6) )
fig11.subplots_adjust(hspace=0.5)
#subplot(3,1,1)
axs11[0].plot( tss, Jleak[:,0], mycolor )
axs11[0].legend(['JRyR$_{tot}$'], fontsize=20 );
#subplot(3,1,2)
axs11[1].plot( tss, Jleak[:,1], mycolor )
axs11[1].legend(['Passive Leak'], fontsize=20 );
#subplot(3,1,3)
axs11[2].plot( tss, Jleak[:,0]-Jleak[:,1], mycolor )
axs11[2].legend(['SR Ca release'], fontsize=20 );
#
#
#%saveas(gcf, 'RyR', 'fig');
#saveas(gcf, [fullPath, 'RyR'], 'png');
fig11.savefig( dataDir + '/RyR.png' )

#
# 
#
#%APD90
#
#load apds_1Hz.txt
#APD90 = apds_1Hz(end,2) % last beat (500th beat)
#
# plt.show()
