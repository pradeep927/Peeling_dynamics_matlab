clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   SIMULATION OF SOFT ADHESION MEDIATED BY MOBILE BINDERS:     %%
%%              2D VARIATIONAL MODELLING                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The 2D system can be understood as a vesicle of length L0 amd aera A0
% subject to a fixed tension T by a loading device (micropipette) with
% mobile binders attached to a symetric vesicle thanks to bound molecules
% through a patch of size s. The contact angle is theta on the patch and
% beta on the micropipette
% The binding kinetics are governed by the two rates kon and koff.
% The force excerted by the pipette is F.
% The diffusion constant is D1 for bonds and D2 for free binders.
% Ntot is known here, the tension is determined by it.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    PLAN    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% 0- Choose type of bonds:       0-> Ideal bonds      1-> Slip bonds
%
% 1- Parameters definition and non dimensionalization
%
% 2- Initialization
%
%       - Geometry -> R,s
%       - Concentrations (equilibrium) -> u0,v0,w0
%       - T,P,F,A0
%       - Discretization (main-b,main-u,mini) + Bsplines
%       - Stiffness matrices -> K,M,A + Bv,Bw + L
%       - Mini-model -> u_g, du_g
%       - Time dicretization parameters
%
% 3- Loop
%
%       - Compute the current value of v for the linearization -> vn
%       - Form the special matrices -> Mnl, Moff and A2
%       - Cu,Cv,Cw
%       - MM, NN
%       - (MM,NN)+L for BC's
%       - Solve the matrix system -> X -> U,V,W
%       - Update s_dot -> update s
%       - Solve the mechanics -> theta, beta
%       - Track the balance of mass
%       - Mini-model -> u_g, du_g
%
% 4- Intersting plots
%
%       - s, total-mass, shape etc...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Toolbox_main')           % Toolbox for the parent programm
addpath('Toolbox_mechanics')      % Toolbox for the mechanics
addpath('Toolbox_bsplines')       % Toolbox for the bsplines calculations
addpath('Toolbox_minimodel')      % Toolbox for the mini model programm


% load('last_eq.mat')
load('Initial_equilibrium_mixed.mat')

c_max=20;
y1=0.9;
g=2;
global Mat_param ini_Val Invariant_mini Invariant_main_b Invariant_main_u
Mat_param.aaa=13;
['Invariant_main_b',num2str(Mat_param.aaa+1),'=Invariant_main_b;'];
eval(['Invariant_main_b',num2str(Mat_param.aaa+1),'=Invariant_main_b;'])
r=[1];
% value1     = 'Ideal bonds=0 , slip bonds=1 or catch bonds=2 ?';
% type       = input(value1);
% options    = optimoptions(@fmincon,'Algorithm','interior-point');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Problem definition and non dimensionalization
% 
%  L0n     = 1;     % Normalization factor for lengths: Length of vesicle
%  C0n     = 1;     % Normalization factor for concentrations: Ntot/L0n
%  F0n     = 1;     % Normalization factor for forces: kb*T*C0n
%  T0n     = 1;     % Normalization factor for times: L0^2/D1
% % 
% T       = 1;
% F       = 1;
%T=255pN/mu_m, kappa=0.1 pN mu m 
 tau     =1*24500;%10*406
%  gamma1  = 1737;   % L*sqrt(T/kappa) 
%  gamma2  = 1751;   % L*sqrt(kc0/T)
%  gamma3  = 2165;   % L/x_beta_s
% gamma4  = 3e3;    % L/x_beta_c
%  gamma5  = 4*2157.7; % L/_gamma
% A       = 1;      % k_off_catch/k_off_slip
% 
%   
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set up parameters necessary to compute all the other ones
% 
 K=2;
%  L0              = L0n/L0n;                  % Total length of the vesicle
%  ji              = 0.9;                      % Ratio between Lm and L0
 koff            = tau/T0n;                  % Unbinding rate (s^-1)
 kon             = K*C0n*koff;               % Binding rate (m/(s.mol))
 kappa           = T*F0n*(L0n/gamma1)^2;     % Bending stiffness (N.m)
 k               = (gamma2/L0n)^2*F0n*T/C0n; % Bonds stiffness (N/m)
%  s_dot           = 0;                        % Speed of the interface (m/s)
 D1              = L0n^2/T0n*0.5;                % Diffusionion constant of bonds (m^2/s)
 D2              = 2*D1;                     % Diffion constant of free binders (m^2/s)
 x_beta_s        = L0n/gamma3;               % Height of the ''slip barrier'' for slip bonds
 x_beta_c        = L0n/gamma4;               % Height of the ''catch barrier'' for slip bonds
 x_gamma         = 1.186645e-04*sqrt(1); %L0n/gamma5;
%  theta0          = asin(1/4);                     % Initial contact angle at the interface
%  beta0           = pi/2;                     % Initial contact angle at top
%  f               = 0*0.4;                      % F/T
% 
% 
%  Mat_param.koff     = 0;
%  Mat_param.kon      = 0;
 Mat_param.D1       = D1;
 Mat_param.D2       = D2;
%  Mat_param.x_beta_s = x_beta_s;
%  Mat_param.x_beta_c = x_beta_c;
%  Mat_param.x_gamma  = x_gamma;
%  Mat_param.A        = A;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Initialization
% 
% 
% % Deduce geometric parameters
% 
%  R       = ji/(sin(theta0)+pi-theta0); % Nondimensional radius (relative to L0)
%  s0      = R*sin(theta0);              % Nondimensional extent of adhesion patch
%  s       = s0;
% % 
%  beta    = beta0;
%  theta   = theta0;
% 
% % Initial uniform concentrations u,v,w chosen at the equilbrium
% gg               = Mat_param.kon/Mat_param.koff;
% 
% v_0              = 1/(2*s*gg)*(-1+sqrt(1+4*gg*s));
% w_0              = v_0;
% u_0              = gg*v_0*v_0;
% 
% % Tension and force, tension, pressure chosen as corresponding to the equiliubrium
%  Mat_param.T      = 0.5*u_0/(1-cos(theta0));   % Tension of the vesicle (N)
%  Mat_param.F      = f*Mat_param.T;             % Applied force on the pipette (N)
%  Mat_param.kappa  = Mat_param.T*F0n*(L0n/gamma1)^2;     % Bending stiffness (N.m)
  Mat_param.k      = (gamma2/L0n)^2*F0n*Mat_param.T/C0n*1;
%  Mat_param.P_num  = Mat_param.T/R;             % Pressure computed thanks to Laplace law
% % % Compute the vesicle properties at zero force at the top: area and total length
%  L_mech           = Length(theta0, beta0, s);  % Total length (without the ''aspired''part)
%  Mat_param.A0     = Area(theta0, beta0, s);    % Area of the vesicle
%  H0               = Height(theta0, beta0, s);  % height of the vesicle until the pipette
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % set up discretization
% ini_Val.p = 2;
% ini_Val.d = 2;
% ini_Val.gaussno = 2;
ini_Val.nconstraints=5;
% 
% % Dicretization of the bond,unbound and mini region and compute the basis
% % functions and their relatives
% 
%  Invariant_main_b.num_el=250;     % Number of elements in the bound region  150
%  Invariant_main_b.num_cp=Invariant_main_b.num_el+2;     % Number of elements in the bound region  152
%  Invariant_main_u.num_el=500;    % Number of elements in the unbound region
  discrete_pre_cal2
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Create stiffness matrices and Lagrange multipliers matrices
% 
% % Initialization of the concentration vector
%  UU    = ones(Invariant_main_b.num_cp,1)*u_0;
%  VV    = ones(Invariant_main_b.num_cp,1)*v_0;
%  WW    = ones(Invariant_main_u.num_cp,1)*w_0;
% kon=0;
% koff=0;
num=size(WW,1);
 WW=max(WW)*ones(num,1);

%%%%%%%%%%%%%%%%%
%%%%%%%%%time parameters
tt=clock;
time1=(5*24-2)*3600;

%%%%%%%%%%%%%%%%%%
type=0;
% 
 X_= [UU; VV; WW];
% 
Mat_param.koff=koff;
Mat_param.kon=kon;
% % % Create the sitiffness matrices for each region
% Bond region
[K_b, M_b, A_b] = create_matrix(Invariant_main_b);
% Unbound region
[K_u, M_u, A_u] = create_matrix(Invariant_main_u);

% Auxilliary matrices
[Bv,Bw,L] = auxiliary_matrices(ini_Val,Invariant_main_b,Invariant_main_u);

% Zero matrices for the total matrix
zeros_b_u = zeros(Invariant_main_b.num_cp,Invariant_main_u.num_cp) ;
zeros_u_b = zeros(Invariant_main_u.num_cp,Invariant_main_b.num_cp) ;
zeros_b_b = zeros(Invariant_main_b.num_cp,Invariant_main_b.num_cp) ;
% 
% % Create a vector for the pressure
% P_                 = ones(1,2*Invariant_main_b.num_el);
% Invariant_mini.P_  = [P_,zeros(1,size(Invariant_mini.gaussp,2)-size(Invariant_main_b.gaussp,2))];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Solve the mini model
% 
  % Slip and catch bonds
    
    % Find the separation distance u_g_ (and its derivate du_g_) minimizing the problem
    [u_g_,du_g_,u_g,u_g_1,du_g,du_g_1,u_g_0,du_g_0]  = sub_scale_slip_catch_complete(Mat_param,Invariant_main_b,UU,theta,s);
 
% plot(1:size(u_g,1),u_g*35000);
    

    % Store the values of the separation
    sep              = [];
    sep              = [sep,u_g];
    dsep             = [];
    dsep             = [dsep,du_g];
    
% u_g_0c=Invariant_main_b.Nsh_b_0(1,:)*ZZ(1:3);
% du_g_0c=Invariant_main_b.Nsh_b_0(2,:)*ZZ(1:3);
% u_g_1c=Invariant_main_b.Nsh_b_0(1,:)*ZZ(150:152);
% du_g_1c=Invariant_main_b.Nsh_b_0(2,:)*ZZ(150:152);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over time
% Parameters for the loop
delta_t =2*4.08163265e-7;           % Time-step
num_t = 150000;               % Number of time-steps


% Initialize vectors for s and shape
s_hist=zeros(num_t,1);    % Initialize the vector storing s
X_hist   = [UU;VV;WW];              % Initialize the vector storing the concentrations
time_hist   = [0];  
u_g_hist = [];
   %%%%%%%%%%%%%%%%555 
X_hist1   = [UU]; 
X_hist2   = [VV]; 
X_hist3   = [WW]; 
angle=zeros(num_t,2);  
total_mass=zeros(num_t,4); 

angle(1,1)=theta; 
angle(1,2)=beta; 

dt_hist=[];

rate_uu=[];
rate_vv=[];
rate_ww=[];
pradeep=0;
div_flux_u=[];
div_flux_v=[];
div_flux_w=[];
div_flux_u_diff=[];
s_hist(1)=s;
% 
F1        = 0.3*Mat_param.T;
num_r     = 10;
F_        = F1*ones(num_t,1);
for i = 1:num_r-1
    F_(i) = (F1-0)*(i)/(num_r-1);
end   
timep=0;

for i = 1:num_t;
      Mat_param.F=F_(i);
  i  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Form the matrices to solve the matrix equation
    
    % Interpolate v for the linearized term
    u_nl= interpolate(Invariant_main_b, UU);
    v_nl= interpolate(Invariant_main_b, VV);
    w_nl= interpolate(Invariant_main_u, WW);
    
    % Compute the special matrices (for ideal bonds A2=zeros and M_off=M_b)
    [M_on, M_off, K11, K12, K21, K22, K3, S] = create_special_matrix_slip_catch_crowd(v_nl,u_nl,w_nl,u_g,du_g,type,c_max);
    
    %Compute the total matrices
    Cu = 1/delta_t*M_b + Mat_param.koff*M_off + 1/s^2*Mat_param.D1*K11 - s_dot/s*A_b +Mat_param.D1/(s^2)*(2*S/((Mat_param.x_gamma)^2));
    Cv = 1/delta_t*M_b +Mat_param.kon*M_on + 1/s^2*Mat_param.D2*(K22-Bv*(1+VV(end)/(c_max-UU(end)-VV(end)))) - s_dot/s*A_b;
    Cw = 1/delta_t*M_u + 1/(1-s)^2*Mat_param.D2*(K3+Bw*c_max/(c_max-WW(1))) - s_dot/(1-s)*A_u;
    
    % Create the system of equations MM*u_n+1 =NN*u_n
    MM = [Cu (-Mat_param.kon*M_on+ 1/s^2*Mat_param.D1*K12) zeros_b_u ;...
        (-Mat_param.koff*M_off+ 1/s^2*Mat_param.D2*(K21-Bv*VV(end)/(c_max-UU(end)-VV(end)))) Cv zeros_b_u;...
        zeros_u_b zeros_u_b Cw];
    
    NN = [1/delta_t*M_b zeros_b_b zeros_b_u; ...
        zeros_b_b 1/(delta_t)*M_b zeros_b_u; ...
        zeros_u_b zeros_u_b 1/delta_t*M_u];
    
    
    
    
    % Mechanical equilibrioum
    L(1,Invariant_main_b.num_cp) = 1;             
    % v(1)=w(0)
    L(2,2*Invariant_main_b.num_cp) = c_max;             
    L(2,2*Invariant_main_b.num_cp+1) = -(c_max-UU(end));  

    
    % D2/s * v'(1) = D2/(L0-s) * w'(0)
    L(3,Invariant_main_b.num_cp-2:Invariant_main_b.num_cp) = -Mat_param.D2/s*Invariant_main_b.Nsh_b_1(2,:)*VV(end)/(c_max-UU(end)-VV(end));
    L(3,2*Invariant_main_b.num_cp-2:2*Invariant_main_b.num_cp) = -Mat_param.D2/s*Invariant_main_b.Nsh_b_1(2,:)*(c_max-UU(end))/(c_max-UU(end)-VV(end));
    L(3,2*Invariant_main_b.num_cp+1:2*Invariant_main_b.num_cp+3) = +(Mat_param.D2)/(1-s)*Invariant_main_u.Nsh_u_0(2,:)*(c_max/(c_max-WW(1)));
       
    T=Mat_param.T;
    y_sol = fsolve(@(x)dirichlet_func(x,c_max,theta,T),y1);
    y1=y_sol;
     if ini_Val.nconstraints == 5
        L(4,1:3)=(Invariant_main_b.Nsh_b_0(1,:)*2*du_g_0+u_g_0/((Mat_param.x_gamma)^2) + Invariant_main_b.Nsh_b_0(2,:)*(1+2*UU(1)/(c_max-UU(1)-VV(1))));
        L(4,Invariant_main_b.num_cp+1:Invariant_main_b.num_cp+3)=(Invariant_main_b.Nsh_b_0(2,:)*2*UU(1)/(c_max-UU(1)-VV(1)));
        L(5,Invariant_main_b.num_cp+1:Invariant_main_b.num_cp+3)=( Invariant_main_b.Nsh_b_0(2,:)*(1+VV(1)/(c_max-UU(1)-VV(1))));
        L(5,1:3)=(Invariant_main_b.Nsh_b_0(2,:)*VV(1)/(c_max-UU(1)-VV(1)));
        
        y=[y_sol; 0; -(WW(1)-VV(end))*s_dot; 0; 0];
        
    else
        y = [y_sol; 0; -(WW(1)-VV(end))*s_dot];
    end
    
    
    
    
    
    % Full matrices with Langrange multipliers
    MM_ = [MM L'; L zeros(ini_Val.nconstraints)];
    bb = [NN*X_ ; y];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solution and update
    XX = MM_\bb;
    X  = XX(1:end-ini_Val.nconstraints);
    X_ = X;
    
    UU = X(1:Invariant_main_b.num_cp);
    
    % Compute the corresponding s_dot
    if type == 0        % Ideal bonds
        du_at_1= Invariant_main_b.Nsh_b_1(2,:)*X(Invariant_main_b.num_cp-2:Invariant_main_b.num_cp);
        dv_at_1= Invariant_main_b.Nsh_b_1(2,:)*X(2*Invariant_main_b.num_cp-2:2*Invariant_main_b.num_cp);
        s_dot = -Mat_param.D1/(s*UU(end))*(UU(end)*2*du_g_1*u_g_1/((Mat_param.x_gamma)^2)+ du_at_1+2*UU(end)*(du_at_1+dv_at_1)/(c_max-UU(end)-VV(end)));
    elseif type == 1    % Slip bonds
       du_at_1= Invariant_main_b.Nsh_b_1(2,:)*X(Invariant_main_b.num_cp-2:Invariant_main_b.num_cp);
        dv_at_1= Invariant_main_b.Nsh_b_1(2,:)*X(2*Invariant_main_b.num_cp-2:2*Invariant_main_b.num_cp);
        s_dot = -Mat_param.D1/(s*UU(end))*(UU(end)*2*du_g_1*u_g_1/((Mat_param.x_gamma)^2)+ du_at_1+2*UU(end)*(du_at_1+dv_at_1)/(c_max-UU(end)-VV(end)));
    else                % Cactch bonds
        du_at_1= Invariant_main_b.Nsh_b_1(2,:)*X(Invariant_main_b.num_cp-2:Invariant_main_b.num_cp);
        f_1    = exp(u_g_1/Mat_param.x_beta_s)+Mat_param.A*exp(-u_g_1/Mat_param.x_beta_c);
        g_1    = exp(u_g_1/Mat_param.x_beta_s)-Mat_param.A*(Mat_param.x_beta_s/Mat_param.x_beta_c)*exp(-u_g_1/Mat_param.x_beta_c);
        g_1    = g_1/f_1;
        s_dot = -Mat_param.D1/(s*UU(end))*(UU(end)*du_g_1*g_1/Mat_param.x_beta_s+ du_at_1);
        
    end
    
    % Update s
    s = s + s_dot * delta_t
    
    % Store the results s, s_dot and the shape
    s_hist(i)=s;
    s_dot_hist(i)=s_dot;
%     X_hist = [X_hist X];
    timep=timep+delta_t;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the mechanics to update theta and beta
    
    % Initial guess for fmincon
    X0(1) = theta;
    X0(2) = beta;
    
    % Solve for vesicle shape
    [Xmec E_val exitflag output lambda] = fmincon(@(X) Energy(X,s,Mat_param.T,Mat_param.F), X0, [], [], [], [], [], [], @(X) constr_eq(X,s,Mat_param.A0),options);
    theta = Xmec(1);
    beta  = Xmec(2);
    
    % Plot the shape
    %plot_conf(theta, beta, s, 'k-')
    Mat_param.P_num = (Mat_param.T/R);
    
    % Store theta and beta
    angle(i,1)=theta;
    angle(i,2)=beta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Track balance of mass
    UU = X(1:Invariant_main_b.num_cp);
    VV = X(Invariant_main_b.num_cp+1:2*Invariant_main_b.num_cp);
    WW = X(2*Invariant_main_b.num_cp+1:end);
    
    
         %%%%%%%%%%%%%%%%%%%%%%%5
       X_hist1 = [X_hist1 UU];
       X_hist2 = [X_hist2 VV];
       X_hist3 = [X_hist3 WW];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Int1 = integrate(Invariant_main_b, UU);
    Int2 = integrate(Invariant_main_b, VV);
    Int3 = integrate(Invariant_main_u, WW);
    total_mass(i,1)  = s*(Int1+Int2) + (1-s)*Int3;   % Total_mass
    total_mass(i,2)  = s*Int1;                       % Number of bonds
    total_mass(i,3)  = s*Int2;                       % Number of free binders in the patch
    total_mass(i,4)  = (1-s)*Int3;                   % Number of free binders in the free part
    total_mass(i,1)/total_mass(1,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remeshing
%      if (s<s_hist(1)*0.3) && (g==1)
%          Mat_param.aaa=5
%          remeshing_continue2
%          g=2;
%          r=[r,i];
%      else
%      end
%      if (s<s_hist(1)*0.15) && (g==2)
%          Mat_param.aaa=2
%          remeshing_continue2
%          g=3;
%          r=[r,i];
%      else
%      end
%      
%      if (s<s_hist(1)*0.002) && (g==3)
%          Mat_param.aaa=1
%          remeshing_continue2
%          g=4;
%          r=[r,i];
%      else
%      end
     
     X_= [UU; VV; WW];

% Create the sitiffness matrices for each region
% Bond region
[K_b, M_b, A_b] = create_matrix(Invariant_main_b);
% Unbound region
[K_u, M_u, A_u] = create_matrix(Invariant_main_u);

% Auxilliary matrices
[Bv,Bw,L] = auxiliary_matrices(ini_Val,Invariant_main_b,Invariant_main_u);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the mini model
    
    if type ~=5    % Slip/catch bonds
        
    [u_g_,du_g_,u_g,u_g_1,du_g,du_g_1,du_g_0]  = sub_scale_slip_catch_complete(Mat_param,Invariant_main_b,UU,theta,s);
        
        % Store the values of the separation
        sep              = [sep,u_g];
        dsep             = [dsep,du_g];
    else            % Ideal bonds
    end
  
  u_g_hist =[u_g_hist,u_g];
   %%%%555%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 



 dt_hist=[dt_hist, delta_t*2450];
 time_hist=[time_hist timep*2450/25];  
    
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   time start
    
    time_elapsed=etime(clock,tt);
    
    
    if time_elapsed>time1 
        break;
    end
  disp('Time elapsd in Hr:');  
 time_elapsed/3600      
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  time end
    

    
    %%%%%%%%%%%%%%%%
   

  if rem(i,2000)==0 || rem(i,500)==0
                save('Res_ideal_cmax20_f03.mat');
  end
     
    if rem(i,50)==0 && s<0.0008
       save('Res_ideal_cmax20_f03_v1.mat');
    end
  

    
    if delta_t < 8.1632650e-06*1.0
    delta_t=1.005*delta_t;
   end

    if delta_t> 8.163265e-6*1.0
       delta_t=8.163265e-6*1.0;
    end
    
   if  s<0.004
     delta_t=0.99*delta_t;
      
   if delta_t < 8.1632650e-07*2
    delta_t=1*8.1632650e-07*2;
    end
    
   end

    
   
    
    disp('time step');
    delta_t*2450
    
    disp('Total time');
    time_hist(i)
    
     if abs(s)>1.0e+4
          save('Res_ideal_cmax20_f03.mat');
        break;
    end
    
end

r=[r,num_t];
