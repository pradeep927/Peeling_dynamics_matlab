function  [dirichlet] = dirichlet_func(x,c_max,theta,T)
global Mat_param



dirichlet=T*(1-cos(theta))+c_max*log((c_max-x)/c_max)+0.5*x;