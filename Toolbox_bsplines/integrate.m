function Integral = integrate(Invariant, U)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate the functions over the domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global ini_Val 
gaussno = ini_Val.gaussno;
Nshape  = Invariant.Nshape;
connect = Invariant.connect;
num_el  = Invariant.num_el;
gweight = Invariant.gweight;


Integral = 0;
for i = 1:num_el
    for k = 1:gaussno
        gp_global = (i-1)*gaussno + k;
        N_Matrix = Nshape(gp_global,:);
        connect_loc = connect(i,:);
        Integral = Integral + N_Matrix * U(connect_loc) *  gweight(gp_global);
    end
end



