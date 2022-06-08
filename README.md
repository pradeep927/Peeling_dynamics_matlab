# Peeling_dynamics_matlab
It contains the MATLAB code for SIMULATION OF SOFT ADHESION MEDIATED BY MOBILE BINDERS
This folder contains several Matlab scripts implementing the computational method described in Supplementary Note 1 to approximate the mathematical model described in the paper. These codes were used to generate all the data in the paper. 

1. diff_F02_rigid_main.m corresponds to diffusion-dominated regime with rigid bonds for F=0.2 (see Fig. 2).
2. diff_F02_compliant_main.m corresponds to diffusion-dominated regime with compliant bonds for F=0.2 (see Fig. 2).
3. reaction_t32_d1_f1_main.m corresponds to reaction-dominated regime (see Fig. 3).
4. mixed_F02_main.m corresponds to mixed regime with ideal bonds for F=0.2, but it can be modified to consider slip bonds (see Fig. 4).
5. crowding_cmax20_F03_main.m corresponds to mixed regime with molecular crowding  (see Fig. 5).

The Matlab Initial*.mat files contain the initial state of the system at the beginning of the simulation.
