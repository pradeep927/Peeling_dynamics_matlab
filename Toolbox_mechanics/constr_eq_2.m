function [c ceq] = constr_eq_2(X, s, A0, H0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint for are conservation in the mechanical model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%c = [];

theta = X(1);
beta = X(2);

c1 = -X(1);
c2 = X(1)-pi/2;
c3 = -X(2);
c4 = X(2)-pi;

ceq1 = Area(theta, beta, s)-A0;
ceq2 = Height(theta, beta, s)-H0;


c   = [c1,c2,c3,c4];
ceq = [ceq1; ceq2];
