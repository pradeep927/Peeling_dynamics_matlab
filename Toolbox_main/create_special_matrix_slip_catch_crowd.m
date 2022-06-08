function [M_on, M_off, K11, K12, K21, K22, K3, S] = create_special_matrix_slip_catch_crowd(v_nl,u_nl,w_nl,u_g,du_g,type,c_max);
global ini_Val  Mat_param Invariant_main_b Invariant_main_u

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
M_on    = zeros(num_cp);    % Create matrix for the linearized term
M_off   = zeros(num_cp);    % Create matrix for the unbinding term 1
K11     = zeros(num_cp);    % Create matrix for the unbinding term 2
K12     = zeros(num_cp);    % Create matrix for the unbinding term 2
K21     = zeros(num_cp);    % Create matrix for the unbinding term 2
K22     = zeros(num_cp);
S       = zeros(num_cp);



% Multiply at every gausspoints
for j = 1:num_el
    M_on_loc  = zeros(size_loc);
    M_off_loc = zeros(size_loc);
    K11_loc   = zeros(size_loc);
    K12_loc   = zeros(size_loc);
    K21_loc   = zeros(size_loc);
    K22_loc   = zeros(size_loc);
    S_loc     = zeros(size_loc);
    
    
    for k = 1:gaussno
        gp_global = (j-1)*gaussno + k;
        N_Matrix = Nshape(gp_global,:);
        B_Matrix = dNsh(gp_global,:);
        connect_loc = connect(j,:);
        M_on_loc = M_on_loc + N_Matrix'*N_Matrix *v_nl(gp_global)*exp(-(u_g(gp_global)/Mat_param.x_gamma)^2)*gweight(gp_global);
        S_loc    = S_loc + B_Matrix'*N_Matrix * du_g(gp_global) * u_g(gp_global)* gweight(gp_global);
        K11_loc  = K11_loc + B_Matrix'*B_Matrix* gweight(gp_global) *(1+ 2*u_nl(gp_global)/(c_max-u_nl(gp_global)-v_nl(gp_global)));
        K12_loc  = K12_loc + B_Matrix'*B_Matrix* gweight(gp_global) *2*u_nl(gp_global)/(c_max-u_nl(gp_global)-v_nl(gp_global));
        K21_loc  = K21_loc + B_Matrix'*B_Matrix* gweight(gp_global) *v_nl(gp_global)/(c_max-u_nl(gp_global)-v_nl(gp_global));
        K22_loc  = K22_loc + B_Matrix'*B_Matrix* gweight(gp_global) *(1+ v_nl(gp_global)/(c_max-u_nl(gp_global)-v_nl(gp_global)));
        
        if type == 1      % Slip bonds
            M_off_loc= M_off_loc + N_Matrix'*N_Matrix * gweight(gp_global)*exp(u_g(gp_global)/Mat_param.x_beta_s);
            
            
        elseif type == 2  % Catch bonds
            
        else              % Ideal bonds
            M_off_loc= M_off_loc + N_Matrix'*N_Matrix * gweight(gp_global);
        end
        
    end
    
    % Assembly
    
    connect_loc = connect(j,:);
    M_on(connect_loc,connect_loc)=  M_on(connect_loc,connect_loc) + M_on_loc;
    M_off(connect_loc,connect_loc)=  M_off(connect_loc,connect_loc) + M_off_loc;
    S(connect_loc,connect_loc)=  S(connect_loc,connect_loc) + S_loc;
    K11(connect_loc,connect_loc)= K11(connect_loc,connect_loc) + K11_loc;
    K12(connect_loc,connect_loc)= K12(connect_loc,connect_loc) + K12_loc;
    K21(connect_loc,connect_loc)= K21(connect_loc,connect_loc) + K21_loc;
    K22(connect_loc,connect_loc)= K22(connect_loc,connect_loc) + K22_loc;
    
   
    
end

Nshape  = Invariant_main_u.Nshape;
dNsh    = Invariant_main_u.dNsh;
connect = Invariant_main_u.connect;
num_cp  = Invariant_main_u.num_cp;
gaussp  = Invariant_main_u.gaussp;
region  = Invariant_main_u.region;
gweight = Invariant_main_u.gweight;
num_el  = Invariant_main_u.num_el;

size_loc = size(connect,2);
K3    = zeros(num_cp);    % Create matrix for the linearized term

for j = 1:num_el
    
    K3_loc   = zeros(size_loc);
    
    
    for k = 1:gaussno
        gp_global = (j-1)*gaussno + k;
        N_Matrix = Nshape(gp_global,:);
        B_Matrix = dNsh(gp_global,:);
        connect_loc = connect(j,:);
        K3_loc  = K3_loc + B_Matrix'*B_Matrix* gweight(gp_global) *(c_max)/(c_max-w_nl(gp_global));
        
        
        
    end
    
    % Assembly
    
    connect_loc = connect(j,:);
    
    K3(connect_loc,connect_loc)= K3(connect_loc,connect_loc) + K3_loc;
    
    
   
    
end
