function [BdNshape BddNshape] = augment_Invariant(Invariant,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Augment Invariant-mini adding extended B-splines vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BdNshape = zeros(size(Invariant.BNshape));
BddNshape = zeros(size(Invariant.BNshape));
for i = 1:length(Invariant.gaussp)
    BdNshape(i,Invariant.spv(i)+[0:p]) =  Invariant.dNsh(i,:);
    BddNshape(i,Invariant.spv(i)+[0:p]) =  Invariant.ddNsh(i,:);
end
