function discrete_pre_cal
global ini_Val Invariant_main_b Invariant_main_u Invariant_mini

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dicretize the two regions of the vesicle  'main_b' and 'main-u' and compute 
% the corresponding shape function and their derivates. Evalueate the
% b-splines at some interesting points (0 and 1).                 
% Store all intersting infos in Invariant_main_b and Invariant_main_u.
% 
% Discretize the extended bond region 'mini' and compute the corresponding
% b-spline and its derivates. Evaluate the bsplines at the extremities of 
% the bonds region (coressponding to 0 and 1/alpha for teh extended domain)
% Store all interesting infos in Invariant_mini
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p       = ini_Val.p;
de      = ini_Val.d;
gaussno = ini_Val.gaussno;

% Dicretization of the bond region
num_el_b = Invariant_main_b.num_el; % Number of elements

U_b = [0:1/num_el_b:1];
U_b = (U_b+10^-13).^(1/10^13);
U_b = U_b-U_b(1);
U_b = U_b/U_b(end);
U_b0= U_b;
h_max_b = max(diff(U_b));
list_el_b = 1:num_el_b;
aux = ones(1,p);
U_b = [aux*U_b(1) U_b aux*U_b(end)];

% Dicretization for the unbound region
num_el_u =  Invariant_main_u.num_el;                            % Number of elements
U_u = [0:1/num_el_u:1];
U_u = (1-U_u+0.001).^(1/100);
U_u = U_u-U_u(1);
U_u = U_u/U_u(end);
h_max_u = max(diff(U_u));
list_el_u = 1:num_el_u;
U_u = [aux*U_u(1) U_u aux*U_u(end)];

% Compute the basis functions 
[Nshape_b,dNsh_b,ddNsh_b,connect_b,spv_b,gweight_b,gaussp_b] = basis_functions(U_b,p,gaussno,de);
[Nshape_u,dNsh_u,ddNsh_u,connect_u,spv_u,gweight_u,gaussp_u] = basis_functions(U_u,p,gaussno,de);

num_cp_b = max(max(connect_b));
num_cp_u = max(max(connect_u));

% Evaluate the derivate for BC D2/s * v'(1) = D2/(L0-s) * w'(0)
[Nsh_b_1,sp]=evaluate_BSp(U_b,p,de,1);   
[Nsh_u_0,sp]=evaluate_BSp(U_u,p,de,0);
[Nsh_b_0,sp]=evaluate_BSp(U_b,p,de,0);   
[Nsh_u_1,sp]=evaluate_BSp(U_u,p,de,1);

% Discretize the mini model region according the the main model

alpha                = 1.4;  % Total length of the mini model region L_mini = alpha*s
U_mini               = [U_b0,1+flip(1-U_b0(1:end-1))];
U_mini(U_mini>alpha) = [];
U_mini(end)          = alpha;
U_mini               = U_mini/alpha;
num_el_m             = size(U_mini,2)-1;
U_mini               = [aux*U_mini(1) U_mini aux*U_mini(end)];
U0                   = zeros(size(U_mini,2)-3,1)';
num_mini             = size(U_mini,2)-3;

[Nshape_m,dNsh_m,ddNsh_m,connect_m,spv_m,gweight_m,gaussp_m] = basis_functions(U_mini,p,gaussno,de);

num_cp_m = max(max(connect_m));

[Nsh_mini_0,sp]=evaluate_BSp(U_mini,p,de,0);
[Nsh_mini_1_alpha,sp] = evaluate_BSp(U_mini,p,de,1/alpha);

Invariant_main_b.U =U_b;
Invariant_main_b.region =1;
Invariant_main_u.region =2;
Invariant_main_b.Nshape = Nshape_b;
Invariant_main_u.Nshape = Nshape_u;
Invariant_main_b.dNsh = dNsh_b;
Invariant_main_u.dNsh = dNsh_u;
Invariant_main_b.ddNsh = ddNsh_b;
Invariant_main_u.ddNsh = ddNsh_u;
Invariant_main_b.spv = spv_b;
Invariant_main_u.spv = spv_u;
Invariant_main_b.gweight = gweight_b;
Invariant_main_u.gweight = gweight_u;
Invariant_main_b.connect=connect_b;
Invariant_main_u.connect=connect_u;
Invariant_main_b.gaussp = gaussp_b;
Invariant_main_u.gaussp = gaussp_u;
Invariant_main_b.Nsh_b_1 = Nsh_b_1;
Invariant_main_b.Nsh_b_0 = Nsh_b_0;
Invariant_main_u.Nsh_u_1 = Nsh_u_1;
Invariant_main_u.Nsh_u_0 = Nsh_u_0;
Invariant_main_b.num_cp = num_cp_b;
Invariant_main_u.num_cp = num_cp_u;

Invariant_mini.num_el = num_el_m;
Invariant_mini.alpha  = alpha;
Invariant_mini.Nshape = Nshape_m;
Invariant_mini.dNsh = dNsh_m;
Invariant_mini.ddNsh = ddNsh_m;
Invariant_mini.spv = spv_m;
Invariant_mini.gweight = gweight_m;
Invariant_mini.connect=connect_m;
Invariant_mini.gaussp = gaussp_m;
Invariant_mini.num_cp = num_cp_m;
Invariant_mini.U_mini   = U_mini';
Invariant_mini.U0       = U0';
Invariant_mini.Nsh_mini_0 = Nsh_mini_0;
Invariant_mini.Nsh_mini_1_alpha =Nsh_mini_1_alpha;
Invariant_mini.num_mini   = num_mini;

% Numerical preliminaries
ini_Val.ncontrol = size(U_mini,2)-3;
num_dofs = ini_Val.ncontrol;

%Compute the BNshape and derivates;
ngp=length(Invariant_mini.spv);
ncp=ini_Val.ncontrol;
spvm=Invariant_mini.spv'*ones(1,ini_Val.p+1)+ones(ngp,1)*(0:ini_Val.p);
BNshape=zeros(ngp,ncp);
for g=1:ngp
    BNshape(g,spvm(g,:))=Invariant_mini.Nshape(g,:);
end
Invariant_mini.BNshape = BNshape;

[BdNshape BddNshape] = augment_Invariant(Invariant_mini,ini_Val.p);
Invariant_mini.BdNshape = BdNshape;
Invariant_mini.BddNshape = BddNshape;


