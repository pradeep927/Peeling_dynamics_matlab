function Aux = auxiliary_calculations_rescale(Invariant,p,C,C1,C2,alpha,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute auxiliary useful terms for the calculation of the energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gweight=Invariant.gweight;
Nshape=Invariant.Nshape;
dNsh=Invariant.dNsh;
ddNsh=Invariant.ddNsh;

Aux.dl = sqrt((alpha*s)^2 + C1.^2);
Aux.curv = C2 ./ (Aux.dl.^3);

Aux.t1 = C2 ./ (Aux.dl.^6);
Aux.t2 = -3 * C1 .* (C2.^2) ./ (Aux.dl.^8);

Aux.t3 = C2 ./ (Aux.dl.^5);
Aux.t4 = -5/2 * C1 .* (C2.^2) ./ (Aux.dl.^7);


