function [C,C1,C2] = interpolate_func(Invariant,x)
global ini_Val 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the vector to have values at gausspoints for u and its
% derivates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = ini_Val.p;
Nshape = Invariant.Nshape;
dNsh = Invariant.dNsh;
ddNsh = Invariant.ddNsh;
spv = Invariant.spv;
ngp_total=length(spv);

C = zeros(ngp_total,1);
C1 = zeros(ngp_total,1);
C2 = zeros(ngp_total,1);
for k=1:(p+1)
    C  = C  + Nshape(:,k).*x((spv(:)+k-1));
    C1 = C1 + dNsh(:,k).*x((spv(:)+k-1));
    C2 = C2 + ddNsh(:,k).*x((spv(:)+k-1));
end
