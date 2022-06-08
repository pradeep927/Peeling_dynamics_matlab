function [ders,sp]=evaluate_BSp(U,p,de,x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the b-splines at the point x between 0 and 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%This function evaluate the values of shape functions and derivatives 
%at a given point
%
m = length(U)-1;
n = m-p-1;
mm = m+1;
pm = p+1;
%
%Computing the span
if x>=U(n+1)-(U(n+1)*1e-9);
    sp = n+1;
else
    b=x>=U;
    b(b==0)=[];
    sp=(length(b));
end
%
%Computing Shape functions
Nshape = zeros(pm,pm);
left = zeros(1,p);
right = zeros(1,p);
Nshape(1,1) = 1;
%
for j = 1:p
    left(j) = (x-U(sp+1-j));
    right(j) = (U(sp+j)-x);
    saved =0;
    r=0;
    while r<j;
        rm = r+1;
        Nshape(j+1,rm) = right(1,r+1)+left(1,j-r);
        temp = Nshape(rm,j)/Nshape(j+1,rm);
        Nshape(rm,j+1) = saved+right(1,r+1)*temp;
        saved = left(j-r)*temp;
        r=rm;
    end
    Nshape(j+1,j+1) = saved;
end
%
%Computation of derivatives
for j = 0:1:p
    jm = j+1;
    ders(1,jm) = Nshape(jm,pm);
end
for r = 0:1:p
    rm = r+1;
    s1 = 1;
    s2 = 2;
    a(1,1) = 1;
    for k = 1:de
        km = k+1;
        d = 0;
        rk = r-k;
        rkm = rk+1;
        pk = p-k;
        pkm = pk+1;
        if r >= k
            a(s2,1) = a(s1,1)/Nshape(pkm+1,rkm);
            d = a(s2,1)*Nshape(rkm,pkm);
        end
        if rk >= -1
            j1 = 1;
        else
            j1 = -rk;
        end
        if (r-1) <= pk
            j2 = k-1;
        else
            j2 = p-r;
        end
        for l = j1:j2
            lm = l+1;
            a(s2,lm) = (a(s1,lm)-a(s1,lm-1))/Nshape(pkm+1,rkm+l);
            d = d+a(s2,lm)*Nshape(rkm+l,pkm);
        end
        if r<= pk
            a(s2,km) = -a(s1,km-1)/Nshape(pkm+1,rm);
            d = d+a(s2,km)*Nshape(rm,pkm);
        end
        ders(km,rm) = d;
        j = s1-1;
        s1 = s2;
        s2 = j+1;
    end
end
%
%Multiplying the derivatives by correction factor
ra = p;
for ka = 1:de
    kam = ka+1;
    for ja = 0:p
        jam = ja+1;
        ders(kam,jam)= ders(kam,jam)*ra;
    end
    ra = ra*(p-ka);
end


