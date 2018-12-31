To compile the program:

 
Using Intel Compiler ICC
icc masterCompute.cpp -o singlecell

./singlecell initial_WTstates.txt stim_param.txt testOutputFolder

__________________________________________________________________
Using GCC
g++ masterCompute.cpp -o singlecell

./singlecell initial_WTstates.txt stim_param.txt testOutputFolder

__________________________________________________________________
plotall.m collects the results from masterCompute.cpp and plots the time course of voltage or currents 



____________________________________________________________________________________
Input settings in stim_param.txt

1 Set Na blocker drug concentration

2 Set Kr blocker ratio

3 Set basic cycle length (BCL)

4 Set how many beats.

5 Set Ligand concentration 

6 Adjust cAMKII activity levels (expression = 'WT', 'OE', or 'KO')
 


____________________________________________________________________________________
initial_WTstates.txt is the initial variables file corresponding to SS_rabbit_varNames.txt (total variables are 206, although some of them are not used in the model)

____________________________________________________________________________________
(1) outputs generated in vm_1Hz.txt 
29	Ca bound to Casqn 
30	Ca SR
31	Na dyad
32	Na sl
33	Na cytosol
34	K 
35	Ca dyad
36	Ca sl
37	Ca cytosol
38	Vm

(2) outputs generated in allresult_1Hz.txt 

I_Ca_store,
I_to_store[0],
I_Na_store,
I_K1_store,
Jserca,
IKs_store,
IKr_store,
Jleak[0],
Jleak[1],
ICFTR,
pars1.Incx 


(3) outputs generated in apds_1Hz.txt

1st column: beat
2nd column: APD90

