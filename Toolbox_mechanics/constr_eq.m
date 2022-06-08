function [c ceq] = constr_eq(X, s, A0)

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
c   = [c1,c2,c3,c4];
ceq = Area(theta, beta, s)-A0;
