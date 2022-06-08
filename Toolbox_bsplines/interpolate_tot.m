function [u_int,du_int,ddu_int] = interpolate_tot(Invariant, U)
global ini_Val 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the vector to have values at gausspoints for u and its
% derivates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gaussno = ini_Val.gaussno;
num_el = Invariant.num_el;
Nshape = Invariant.Nshape;
dNsh   = Invariant.dNsh;
ddNsh  = Invariant.ddNsh;
connect= Invariant.connect;



for i = 1:num_el
    for k = 1:gaussno
        gp_global = (i-1)*gaussno + k;
        N_Matrix = Nshape(gp_global,:);
        B_Matrix = dNsh(gp_global,:);
        D_Matrix = ddNsh(gp_global,:);
        connect_loc = connect(i,:)';
        u_int(gp_global)   =  N_Matrix * U(connect_loc);
        du_int(gp_global)  =  B_Matrix * U(connect_loc);
        ddu_int(gp_global) =  D_Matrix * U(connect_loc);
    end
end



