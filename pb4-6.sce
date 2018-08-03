//Pb-4-6 ISBN:978-1-138-19680-3
clc;clear;mode(-1)
exec('PMM.sci');
//IM7 fiber
EA = 276000//MPa, Table 2.1 ISBN:978-1-138-19680-3
ET = EA
vA = 0.2//Table 2.1 ISBN:978-1-138-19680-3
vT = vA
Vf = 0.591
GA = EA/2/(1+vA); disp(GA)
//8552 Epoxy Table 2.13 ISBN:978-1-138-19680-3
Vm = 1-Vf; disp(Vm)
Em = 4667//MPa
vm = 0.38
Gm = Em/2/(1+vm); disp(Gm)
// composite
//E2 = 1/(Vf/EA + Vm/Em); disp(E2) // rule of mixtures->not accurate
[E1,E2,G12,v12,v23,G23] = PMM(EA,ET,GA,vA,vT,Em,vm,Vf); //best predictions
disp([E2,v23,G23])
disp(E2/2/(1+v23))//check transverse shear value 

