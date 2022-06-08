function [M_nl1, M_nl2, M_off, A2, K2] = create_special_matrix_slip_catch_nl(Invariant,v_nl,u_nl,du_nl,u_g,du_g,c_max,type);
global ini_Val  Mat_param Invariant_main_b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the three special matrices corresponding to the linearized
% term, the off rate and the new term coming from non uniformity of the
% chemical potential over the patch
% For ideal bonds -> Mnl, Moff=M_b, A2= zeros
% For slip bonds  -> Mnl, Moff, A2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gaussno = ini_Val.gaussno;

x_beta_s  = Mat_param.x_beta_s;

Nshape  = Invariant_main_b.Nshape;
dNsh    = Invariant_main_b.dNsh;
connect = Invariant_main_b.connect;
num_cp  = Invariant_main_b.num_cp;
gaussp  = Invariant_main_b.gaussp;
region  = Invariant_main_b.region;
gweight = Invariant_main_b.gweight;
num_el  = Invariant_main_b.num_el;

size_loc = size(connect,2);
M_nl1    = zeros(num_cp);    % Create matrix for the linearized term for vn+1
M_nl2    = zeros(num_cp);     % Create matrix for the linarized term for un+1'
M_off    = zeros(num_cp);    % Create matrix for the unbinding term (to change cause of new factor)
A2       = zeros(num_cp);    % Create matrix for the chemical potential term
K2       = zeros(num_cp);     % Create matrix for the new stiffness matrix K (to change cause of new factor)



% Multiply at every gausspoints
for j = 1:num_el
    M_loc2= zeros(size_loc);
    M_loc3= zeros(size_loc);
    A_loc3= zeros(size_loc);
    A_loc2= zeros(size_loc);
    K_loc2= zeros(size_loc);
    
    for k = 1:gaussno
        gp_global = (j-1)*gaussno + k;
        N_Matrix = Nshape(gp_global,:);
        B_Matrix = dNsh(gp_global,:);
        connect_loc = connect(j,:);
        M_loc2= M_loc2 + N_Matrix'*N_Matrix *v_nl(gp_global)* gweight(gp_global);
        A_loc3= A_loc3 + B_Matrix'*N_Matrix * c_max*du_nl(gp_global)/((c_max-u_nl(gp_global))^2)* gweight(gp_global);
        K_loc2= K_loc2 + B_Matrix'*B_Matrix* gweight(gp_global) * c_max/(c_max-u_nl(gp_global));
        
        if type == 1      % Slip bonds
            M_loc3= M_loc3 + N_Matrix'*N_Matrix * gweight(gp_global)*exp(u_g(gp_global)/Mat_param.x_beta_s)*(c_max/(c_max-u_nl(gp_global)));
            A_loc2= A_loc2 + B_Matrix'*N_Matrix * du_g(gp_global) * gweight(gp_global);
        elseif type == 2  % Catch bonds
            f_c = exp(u_g(gp_global)/Mat_param.x_beta_s)+Mat_param.A*exp(-u_g(gp_global)/Mat_param.x_beta_c);
            g_c = exp(u_g(gp_global)/Mat_param.x_beta_s)-Mat_param.A*(Mat_param.x_beta_s/Mat_param.x_beta_c)*exp(-u_g(gp_global)/Mat_param.x_beta_c);
            g_c = g_c/f_c;
            M_loc3= M_loc3 + N_Matrix'*N_Matrix * gweight(gp_global)*f_c*(c_max/(c_max-u_nl(gp_global)));
            A_loc2= A_loc2 + B_Matrix'*N_Matrix * du_g(gp_global)*g_c*gweight(gp_global);  
        else              % Ideal bonds
            M_loc3= M_loc3 + N_Matrix'*N_Matrix * gweight(gp_global)*(c_max/(c_max-u_nl(gp_global)));
        end
        
    end
    
    % Assembly
    
    connect_loc = connect(j,:);
    M_nl1(connect_loc,connect_loc) =  M_nl1(connect_loc,connect_loc) + M_loc2;
    M_nl2(connect_loc,connect_loc) =  M_nl2(connect_loc,connect_loc) + A_loc3;
    K2(connect_loc,connect_loc)    =  K2(connect_loc,connect_loc) + K_loc2;
    
    if type == 1     % Slip bonds
        M_off(connect_loc,connect_loc)=  M_off(connect_loc,connect_loc) + M_loc3;
        A2(connect_loc,connect_loc)=  A2(connect_loc,connect_loc) + A_loc2;
    elseif type == 2     % Catch bonds
        M_off(connect_loc,connect_loc)=  M_off(connect_loc,connect_loc) + M_loc3;
        A2(connect_loc,connect_loc)=  A2(connect_loc,connect_loc) + A_loc2;
    else             % Ideal bonds
        M_off(connect_loc,connect_loc)=  M_off(connect_loc,connect_loc) + M_loc3;
    end
    
end

