clear all

global Mat_param;
y1    =1.5741;
theta =0.2527;
T     = 24.7861;
T_theta=T*(1-cos(theta));


c_max_1=20000;
c_max_2=20;
c_max_3=5;

c_max=1000;

n_pts=1000;
x_max=5;
f_1=zeros(n_pts,1);
f_2=zeros(n_pts,1);
f_3=zeros(n_pts,1);
xx=zeros(n_pts,1);
for i=1:n_pts;
x=x_max*i/n_pts;
xx(i,1)=x;
f_1(i,1)=T_theta+c_max_1*log((c_max_1-x)/c_max_1)+0.5*x;
f_2(i,1)=T_theta+c_max_2*log((c_max_2-x)/c_max_2)+0.5*x;
f_3(i,1)=T_theta+c_max_3*log((c_max_3-x)/c_max_3)+0.5*x;
end
x=0;
plot(xx,f_1)
hold on
plot(xx,f_2)
plot(xx,f_3)

syms z
eqn = dirichlet_func(z,c_max,theta,T);
y = fsolve(@(x)dirichlet_func(x,c_max,theta,T),y1);
