function [u_g]  = sub_scale_visualize(Mat,Invariant,UU,theta,s);
global   ini_Val Invariant_mini Mat_param Invariant_main_b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mini model: Minimization of the energy  energy_mini_rescale to compute
% the separation distance between the vesicles. This distance is first
% computed on an extended domain (length alpha*s) thats is rescale to be
% in [0,1]. Finally we cut the solution to be on [0,s]
%
% C_b is the extended value of u on [0,alpha] (0 added to fit the size)
%
% We also compute the values of the derivates at the extremes points 0 and
% 1/alpha (corresponding to 1 in the patch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definition of parameters
F_mini = sin(theta)*Mat_param.T;
alpha  = Invariant_mini.alpha;

% Compute the value of concentration of bonds u_f
u_f                = interpolate(Invariant_main_b, UU); 

% Fit the vector u_f to the size of the mesh of mini_model (add zeros)
u_f                = [u_f,zeros(1,size(Invariant_mini.gaussp,2)-size(Invariant_main_b.gaussp,2))];

% Make the vector for the pressure
Mat_param.P        = 0* Mat_param.P_num*Invariant_mini.P_;

% Set the concentration distribution
Mat_param.C_b= Mat_param.k*u_f';

% Minimization
options_mini = optimoptions('fminunc','GradObj','on','Hessian','off',...
    'DerivativeCheck','off','Display','notify-detailed','Algorithm','quasi-newton','MaxIter',5000,'TolFun', 1e-14, 'TolX', 1e-14);

[U_,Fval,eFlag,output,gradE_1]      = fminunc(@(U)Energy_mini_rescale2(U,alpha,F_mini,s),Invariant_mini.U0,options_mini);

Invariant_mini.U0                   = U_;            
[Cs_,Cs1_,Cs2_]                     = interpolate_tot(Invariant_mini, U_);
u_g_=Cs_';


% Cut to obtain the vectors on the good interval
u_g                = u_g_([1:size(Invariant_main_b.gaussp,2)]);





