function A = Area(theta, beta, s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the area of the vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


R = s/(sin(theta)+cos(beta));

A = R^2 * (3*pi/2-theta-beta+cos(theta)*sin(theta)+cos(beta)*sin(beta)+2*cos(beta)*cos(theta));
A = A/2;