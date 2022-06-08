function [u_g_,du_g_,u_g,u_g_1,du_g,du_g_1,u_g_0,du_g_0]  = sub_scale_slip_catch_rigid(Mat,Invariant,UU,theta,s);
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


% Cut to obtain the vectors on the good interval
u_g                = 0;
du_g               = 0;
du_g_1             = 0;
du_g_0             = 0;




