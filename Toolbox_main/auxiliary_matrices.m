function [Bv,Bw,L] = auxiliary_matrices(ini_Val,Invariant_main_b,Invariant_main_u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the three auxiliary matrices corresponding to the terms appearing
% during the integration by parts Bv and BW and the constant terms of the Lagrange
% multiplier matrix L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nconstraints = ini_Val.nconstraints;
p            = ini_Val.p;

% Boundary terms apprearing with integration by parts
% gives vn+1'(1)
Bv = zeros(Invariant_main_b.num_cp);
Bv(end,Invariant_main_b.num_cp-ini_Val.p:Invariant_main_b.num_cp) = Invariant_main_b.Nsh_b_1(2,:); 

% gives wn+1'(0)
Bw = zeros(Invariant_main_u.num_cp);                      
Bw(1,1:1+ini_Val.p) = Invariant_main_u.Nsh_u_0(2,:);

% Lagrange multiplier 

L=zeros(nconstraints,2*Invariant_main_b.num_cp+Invariant_main_u.num_cp);
% Boundary/jump conditions with Lagrange multipliers
% u(1)=2T(1-cos(theta))  
L(1,Invariant_main_b.num_cp) = 1;             
% v(1)=w(0)
L(2,2*Invariant_main_b.num_cp) = 1;             
L(2,2*Invariant_main_b.num_cp+1) = -1;  
