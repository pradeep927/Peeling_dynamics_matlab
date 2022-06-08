function [Nshape,dNsh,ddNsh,connect,spv,gweight,gaussp] = basis_functions(U,p,gaussno,de)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the basis functions knowing the mnesh and p,d and the number of
% gauss nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a = U';
a(1:p)=[];
a(end-p+1:end)=[];
b=a;
a(length(a)) = [];
b(1) = [];
[gpoint,gweight_local]=gaquad(gaussno);
h = ((b-a)/2);

Nshape = zeros(length(a)*gaussno,(p+1));
dNsh = zeros(length(a)*gaussno,(p+1));
ddNsh = zeros(length(a)*gaussno,(p+1));
gweight=zeros(length(a)*gaussno,1);
connect=zeros(length(a),p+1);
iglobal=0;
for i=1:length(a)
    for k=1:gaussno
        iglobal=iglobal+1;
        gweight(iglobal)=gweight_local(k)*h(i);
        gaussp(iglobal) =((b(i)-a(i))/2)*gpoint(k)+((a(i)+b(i))/2);
        [Nsh,sp]=evaluate_BSp(U,p,de,gaussp(iglobal));
        Nshape(iglobal,:)=Nsh(1,:);
        dNsh(iglobal,:)=Nsh(2,:);
        ddNsh(iglobal,:)=Nsh(3,:);
        spv(iglobal)=sp-p;
    end
    connect(i,:)=[i:i+p];
end


