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
load('Initial_equilibrium_mixed.mat')
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
 type=1;
 L0n     = 1;     % Normalization factor for lengths: Length of vesicle
 C0n     = 1;     % Normalization factor for concentrations: Ntot/L0n
 F0n     = 1;     % Normalization factor for forces: kb*T*C0n
 T0n     = 1;     % Normalization factor for times: L0^2/D1
% 
 T       = 1;
 F       = 1;
  K       = 2;
 tau     = 24500;
 gamma1  = 1737;   % L*sqrt(T/kappa) 
 gamma2  = 1751;   % L*sqrt(k/T)
 gamma3  = 2165;   % L/x_beta_s
 gamma4  = 3e3;    % L/x_beta_c
 gamma5  = 4*2165; % L/_gamma
 A       = 1;      % k_off_catch/k_off_slip
% 
%   
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set up parameters necessary to compute all the other ones
% 
 L0              = L0n/L0n;                  % Total length of the vesicle
 ji              = 0.9;                      % Ratio between Lm and L0
 koff            = tau/T0n;                  % Unbinding rate (s^-1)
 kon             = K*C0n*koff;               % Binding rate (m/(s.mol))
 kappa           = T*F0n*(L0n/gamma1)^2;     % Bending stiffness (N.m)
 k               = (gamma2/L0n)^2*F0n*T/C0n; % Bonds stiffness (N/m)  (gamma2/L0n)^2*F0n*Mat_param.T/C0n;
 s_dot           = 0;                        % Speed of the interface (m/s)
 D1              = 0.5*L0n^2/T0n*0.001;                % Diffusionion constant of bonds (m^2/s)
 D2              = 2*D1;                     % Diffion constant of free binders (m^2/s)
 x_beta_s        = L0n/gamma3*1;               % Height of the ''slip barrier'' for slip bonds
 x_beta_c        = L0n/gamma4;               % Height of the ''catch barrier'' for slip bonds
 x_gamma         = 1.186646e-04;
 
 theta0          = 32/180*pi;                     % Initial contact angle at the interface
  beta0           = pi/2;                     % Initial contact angle at top
% k;

 % f               = 0*0.4;                      % F/T
% 
% 
 Mat_param.koff     =koff;
  Mat_param.kon      =kon;
 Mat_param.D1       = D1;
Mat_param.D2       = D2;
 Mat_param.x_beta_s = x_beta_s;
 Mat_param.x_beta_c = x_beta_c;
 x_gamma=sqrt(1)*x_gamma; %%%%%%%%%%%%%%%%%%%%%%%%%%% important
 Mat_param.x_gamma  = x_gamma;
% Mat_param.A        = A;


% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Initialization
% 
% 
% % Deduce geometric parameters
%%%%%%%%%%%5
length=6;
L0=length/35;
% s=4/35;
s0=s;
%%%%%%%%%%%%%%%

theta_d=32;%in degree
theta=theta_d*pi/180;
H0=(L0-s)*sin(theta);
H=H0;
angle(1,1)=theta; 
angle(1,2)=H; 
%%%%%%%%%%%
theta0=theta;
s_hist(1)=s;
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
%  Mat_param.k      = (gamma2/L0n)^2*F0n*Mat_param.T/C0n;
%  Mat_param.P_num  = Mat_param.T/R;             % Pressure computed thanks to Laplace law
% % Compute the vesicle properties at zero force at the top: area and total length
%  L_mech           = Length(theta0, beta0, s);  % Total length (without the ''aspired''part)
%  Mat_param.A0     = Area(theta0, beta0, s);    % Area of the vesicle
%  H0               = Height(theta0, beta0, s);  % height of the vesicle until the pipette
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % set up discretization
ini_Val.p = 2;
ini_Val.d = 2;
ini_Val.gaussno = 2;
ini_Val.nconstraints=4;
% 
% % Dicretization of the bond,unbound and mini region and compute the basis
% % functions and their relatives
% 
 Invariant_main_b.num_el=350;     % Number of elements in the bound region

 Invariant_main_b.num_cp=Invariant_main_b.num_el+2;
 Invariant_main_u.num_el=700;    % Number of elements in the unbound region
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
%  type=1;
%  save('slipk01.mat');
num=size(WW,1);
 WW=max(WW)*ones(num,1);

%%%%%%%%%%%%%%%%%
%%%%%%%%%time parameters
tt=clock;
time1=(5*24-2)*3600;

%%%%%%%%%%%%%%%%%%
type=1;
% 
u_0= Mat_param.T *2* (1-cos(theta));
UU(end)=u_0;

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
P_                 = ones(1,2*Invariant_main_b.num_el);
Invariant_mini.P_  = [P_,zeros(1,size(Invariant_mini.gaussp,2)-size(Invariant_main_b.gaussp,2))];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Solve the mini model
% 
  % Slip and catch bonds
    
    % Find the separation distance u_g_ (and its derivate du_g_) minimizing the problem
    [u_g_,du_g_,u_g,u_g_1,du_g,du_g_1,u_g_0,du_g_0]  = sub_scale_slip_catch_complete(Mat_param,Invariant_main_b,UU,theta,s);
 

    

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
delta_t =1.0*8.163265e-7;           % Time-step
num_t = 35000;               % Number of time-steps


% Initialize vectors for s and shape
s_hist=zeros(num_t+5,1);    % Initialize the vector storing s


X_hist   = [UU;VV;WW];              % Initialize the vector storing the concentrations
X_hist1   = [UU]; 
X_hist2   = [VV]; 
X_hist3   = [WW]; 

time_hist   = [0];  
u_g_hist = [];
  
angle=zeros(num_t+5,2);  
total_mass=zeros(num_t+5,4); 



dt_hist=[];
%%%%%%%%%%%%
% theta_d=32;%in degree
% theta=theta_d*pi/180;
% H0=(L0-s)*sin(theta);
% H=H0;
% angle(1,1)=theta; 
% angle(1,2)=H; 
% %%%%%%%%%%%
% 
% s_hist(1)=s;
% % 
F1        = 0.3*Mat_param.T;
num_r     = 10;
F_        = F1*ones(num_t,1);
for i = 1:num_r-1
    F_(i) = (F1-0)*(i)/(num_r-1);
end   
timep=0;

for i = 1:num_t;
     Mat_param.F=F_(i);
%      Mat_param.T=Mat_param.F *sin(theta);
  i  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Form the matrices to solve the matrix equation
    
    % Interpolate v for the linearized term
    v_nl= interpolate(Invariant_main_b, VV);
    
    % Compute the special matrices (for ideal bonds A2=zeros and M_off=M_b)
    [M_nl, M_off, A2] = create_special_matrix(Invariant_main_b,v_nl,u_g,du_g,type);
    
    %Compute the total matrices
    Cu = 1/delta_t*M_b + Mat_param.koff*M_off + 1/s^2*Mat_param.D1*K_b - s_dot/s*A_b +Mat_param.D1/(s^2)*(2*A2/((Mat_param.x_gamma)^2));
    Cv = 1/delta_t*M_b +Mat_param.kon*M_nl + 1/s^2*Mat_param.D2*(K_b-Bv) - s_dot/s*A_b;
    Cw = 1/delta_t*M_u + 1/(L0-s)^2*Mat_param.D2*(K_u+Bw) - s_dot/(L0-s)*A_u;
    
    % Create the system of equations MM*u_n+1 =NN*u_n
    MM = [Cu -Mat_param.kon*M_nl zeros_b_u ;...
        -Mat_param.koff*M_off Cv zeros_b_u;...
        zeros_u_b zeros_u_b Cw];
    
    NN = [1/delta_t*M_b zeros_b_b zeros_b_u; ...
        zeros_b_b 1/(delta_t)*M_b zeros_b_u; ...
        zeros_u_b zeros_u_b 1/delta_t*M_u];
    
    % D2/s * v'(1) = D2/(L0-s) * w'(0)
    L(3,2*Invariant_main_b.num_cp-2:2*Invariant_main_b.num_cp) = Mat_param.D2/s*Invariant_main_b.Nsh_b_1(2,:);
    L(3,2*Invariant_main_b.num_cp+1:2*Invariant_main_b.num_cp+3) = -(Mat_param.D2)/(L0-s)*Invariant_main_u.Nsh_u_0(2,:);
    
    if ini_Val.nconstraints == 4
        L(4,1:3)=(Invariant_main_b.Nsh_b_0(1,:)*2*du_g_0*u_g_0/((Mat_param.x_gamma)^2) + Invariant_main_b.Nsh_b_0(2,:));
%       L(4,1:3)=(Invariant_main_b.Nsh_b_0(1,:)*ddelta_0/Mat_param.x_beta + Invariant_main_b.Nsh_b_0(2,:));
        y = [2*Mat_param.T*(1-cos(theta)); 0; 0; 0];
        
    else
        y = [2*Mat_param.T*(1-cos(theta)); 0; 0];
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
        s_dot = -Mat_param.D1/(s*UU(end))*(UU(end)*2*du_g_1*u_g_1/((Mat_param.x_gamma)^2)+ du_at_1);
    elseif type == 1    % Slip bonds
        du_at_1= Invariant_main_b.Nsh_b_1(2,:)*X(Invariant_main_b.num_cp-2:Invariant_main_b.num_cp);
        s_dot = -Mat_param.D1/(s*UU(end))*(UU(end)*2*du_g_1*u_g_1/((Mat_param.x_gamma)^2)+ du_at_1);
    else                % Cactch bonds
        du_at_1= Invariant_main_b.Nsh_b_1(2,:)*X(Invariant_main_b.num_cp-2:Invariant_main_b.num_cp);
        f_1    = exp(u_g_1/Mat_param.x_beta_s)+Mat_param.A*exp(-u_g_1/Mat_param.x_beta_c);
        g_1    = exp(u_g_1/Mat_param.x_beta_s)-Mat_param.A*(Mat_param.x_beta_s/Mat_param.x_beta_c)*exp(-u_g_1/Mat_param.x_beta_c);
        g_1    = g_1/f_1;
        s_dot = -Mat_param.D1/(s*UU(end))*(UU(end)*du_g_1*g_1/Mat_param.x_beta_s+ du_at_1);
        
    end
    
    disp('patch length');
    % Update s
    s = s + s_dot * delta_t
    
    % Store the results s, s_dot and the shape
    s_hist(i+1)=s;
    s_dot_hist(i)=s_dot;
%     X_hist = [X_hist X];
    timep=timep+delta_t;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the mechanics to update theta and beta
    
    % Initial guess for fmincon
%     X0(1) = theta;
%     X0(2) = beta;
    
    % Solve for vesicle shape
%     [Xmec E_val exitflag output lambda] = fmincon(@(H) Energy(H,s,Mat_param.T,Mat_param.F), H0, [], [], [], [], [], [], @(H) constr_eq(H,s,Mat_param.L0),options);
%     theta = Xmec(1);
%     beta  = Xmec(2);
    H=(L0-s)*sin(theta);
    % Plot the shape
    %plot_conf(theta, beta, s, 'k-')
    Mat_param.P_num = (Mat_param.T/R);
    
    % Store theta and beta
    angle(i+1,1)=theta;
    angle(i+1,2)=H;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Track balance of mass
    UU = X(1:Invariant_main_b.num_cp);
    VV = X(Invariant_main_b.num_cp+1:2*Invariant_main_b.num_cp);
    WW = X(2*Invariant_main_b.num_cp+1:end);
    
          %%%%%%%%%%%%%%%%%%%%%%%5
       X_hist1 = [X_hist1 UU];
       X_hist2 = [X_hist2 VV];
       X_hist3 = [X_hist3 WW];
    
    
    Int1 = integrate(Invariant_main_b, UU);
    Int2 = integrate(Invariant_main_b, VV);
    Int3 = integrate(Invariant_main_u, WW);
    
    
    total_mass(i,1)  = s*(Int1+Int2) + (L0-s)*Int3;   % Total_mass
    total_mass(i,2)  = s*Int1;                       % Number of bonds
    total_mass(i,3)  = s*Int2;                       % Number of free binders in the patch
    total_mass(i,4)  = (L0-s)*Int3;                   % Number of free binders in the free part
    
    disp('Total Mass');
    total_mass(i,1)/total_mass(1,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
%     if abs(s)<0.0008
%         pradeep=pradeep+1;
%         if pradeep>4000
%             break;
%         end
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve the mini model
    if type ~=5  
%    if type ==1 || type==2    % Slip/catch bonds
        
    [u_g_,du_g_,u_g,u_g_1,du_g,du_g_1,u_g_0,du_g_0]  = sub_scale_slip_catch_complete(Mat_param,Invariant_main_b,UU,theta,s);
     
   
%     [u_g151,du_g151,ddu_g,U11, V11, W11, d_U,  d_V, d_W, dd_U, dd_V, dd_W]  = sub_scale_slip_catch_927(Mat_param,Invariant_main_b,UU,VV,WW, theta,s);
    
    
    
    
    
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
        save('Res_reac_t32_d1_f1_len.mat'); 
        break;
    end
  disp('Time elapsd in Hr:');  
 time_elapsed/3600      
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  time end
  disp('length');
  (s+H/sin(theta))*35
    
%     save('Res_reac_06_d2_theta32.mat');
%     if abs(s)>s0
%     r=[r,num_t];
%     save('Res_ideal_mixed_03_m1.mat');
%     elseif s==NaN
%     r=[r,num_t];
%     save('Res_ideal_mixed_03_m1.mat');
%     else
%     end
    
    %%%%%%%%%%%%%%%%
   

  if rem(i,1000000)==0 || rem(i,500000)==0
      save(sprintf('Res_reac_initial_t14_%d.mat',i));
  end
     
    
   if rem(i,100)==0 
      save('Res_reac_t32_d1_f1_len.mat');
  end

 if rem(i,20)==0 && s<0.0005
      save('Res_reac_t32_d1_f1_len.mat');
  end


     

    
    if delta_t < 8.1632650e-07*1.0
    delta_t=1.001*delta_t;
   end

    if delta_t> 8.163265e-7*1.0
        delta_t=8.163265e-7*1.0;
    end
    
    if  s/s0<0.3
       delta_t=0.999*delta_t;
       
    if delta_t < 8.1632650e-07*1.0
     delta_t=1.0*8.1632650e-07;
    end
     
    end

     disp('Height');
     H
    disp('Angle');
    theta
    
    disp('time step');
    delta_t*2450
    
    disp('Total time');
    time_hist(i)
    
     if abs(s)>1.0e+4
      save('Res_reac_t32_d1_f1_len.mat');
        break;
    end
    
end

r=[r,num_t];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interesting plots

% Plot the evolution of total_mass

