function  [Ener,grad_E] = Energy_mini_rescale2(U,alpha,F_mini,s)
global Invariant_mini ini_Val Mat_param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the total energy of the system and its corresponding gradients
% for minimization in the parent file sub-scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[C,C1,C2]=interpolate_func(Invariant_mini, U);

% Compute datas for the next calculations dl,t1t2 etc..
Aux= auxiliary_calculations_rescale2(Invariant_mini,ini_Val.p,C,C1,C2,Invariant_mini.alpha,s);


%% Tension
% Energy
E_T = Mat_param.T*sum(Aux.dl .* Invariant_mini.gweight);
% Gradient
F_T = Mat_param.T *(C1./Aux.dl.* Invariant_mini.gweight)' * Invariant_mini.BdNshape;

%% Bending

E_B = 0.5*(alpha*s)^3*Mat_param.kappa* sum(Aux.curv.^2.* Invariant_mini.gweight);
% Gradient
F_B = (Aux.t1.*Invariant_mini.gweight)'* Invariant_mini.BddNshape + ...
      (Aux.t2.*Invariant_mini.gweight)'* Invariant_mini.BdNshape;
F_B = (alpha*s)^3*F_B * Mat_param.kappa;

%% Bonds
% Energy
E_bonds = 0.5*(alpha*s)*sum(Mat_param.C_b.*C.^2.*Invariant_mini.gweight);
% Gradient
F_bonds = (alpha*s)*(Mat_param.C_b.*C.*Invariant_mini.gweight)'* Invariant_mini.BNshape;

%% Pressure
% % Energy
% E_pressure = (alpha*s)*sum(Mat_param.P'.*C.* Invariant_mini.gweight);
% % Gradient
% F_pressure = (alpha*s)*(Mat_param.P'.*Invariant_mini.gweight)'*Invariant_mini.BNshape;

%% Force at the end
% Energy
E_force = -U(end)*F_mini;
% Gradient
F_force = zeros(size(F_T));
F_force(end) =-F_mini;

%% Total energie and gradient
% Gradient
Ener = E_T + E_B + E_bonds  + E_force;
% Gradient
grad_E = F_T + F_B + F_bonds  + F_force;


grad_E = grad_E';
