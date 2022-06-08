function [K, M, A] = create_matrix(Invariant_main);
global ini_Val 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the three matrices K,M and A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gaussno = ini_Val.gaussno;

Nshape  = Invariant_main.Nshape;
dNsh    = Invariant_main.dNsh;
connect = Invariant_main.connect;
num_cp  = Invariant_main.num_cp;
gaussp  = Invariant_main.gaussp;
region  = Invariant_main.region;
gweight = Invariant_main.gweight;
num_el  = Invariant_main.num_el;


%Build stiffness matrix
size_loc = size(connect,2);
M = zeros(num_cp);
K = zeros(num_cp);
A = zeros(num_cp);

for i = 1:num_el
    M_loc = zeros(size_loc);
    K_loc = zeros(size_loc);
    A_loc = zeros(size_loc);
    for k = 1:gaussno
        gp_global = (i-1)*gaussno + k;
        N_Matrix = Nshape(gp_global,:);
        B_Matrix = dNsh(gp_global,:);
        M_loc = M_loc  + N_Matrix'*N_Matrix * gweight(gp_global);
        K_loc = K_loc  + B_Matrix'*B_Matrix * gweight(gp_global);
        if region == 1
            A_loc = A_loc  + N_Matrix'*B_Matrix * gweight(gp_global)...
                *  gaussp(gp_global);
        else
            A_loc = A_loc  + N_Matrix'*B_Matrix * gweight(gp_global)...
                * (1-gaussp(gp_global));
        end
    end
    connect_loc = connect(i,:);
    M(connect_loc,connect_loc) =  M(connect_loc,connect_loc) + M_loc;
    K(connect_loc,connect_loc) =  K(connect_loc,connect_loc) + K_loc;
    A(connect_loc,connect_loc) =  A(connect_loc,connect_loc) + A_loc;
end

