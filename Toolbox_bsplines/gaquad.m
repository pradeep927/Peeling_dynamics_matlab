function [gpoint,gweight]=gaquad(gaussno);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to assign the gauss points and weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if gaussno == 2;
    gweight = [1 1];
    gpoint = [(-1/sqrt(3)) (1/sqrt(3))];
end
if gaussno == 3;
    gweight = [(5/9) (8/9) (5/9)];
    gpoint = [(-sqrt(15)/5) 0 (sqrt(15)/5)];
end
if gaussno == 4;
    gweight = [((18-sqrt(30))/36) ((18+sqrt(30))/36) ((18+sqrt(30))/36) ((18-sqrt(30))/36)];
    gpoint = [-(sqrt(525+sqrt(30)*70)/35) -(sqrt(525-sqrt(30)*70)/35) (sqrt(525-sqrt(30)*70)/35) (sqrt(525+sqrt(30)*70)/35)];
end
if gaussno == 6;
    [x, w] = lgwt(gaussno,-1,1);
    gpoint = x';
    gweight = w';
end
if gaussno > 10;
    gweight = ones(gaussno);
    gpoint = [-1:(2/(gaussno-1)):1] * 0.999999;
end
