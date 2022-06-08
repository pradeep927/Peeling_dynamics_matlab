function E = Energy(X, s, T, F)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the total energy of the vesicle. Will be minimize uin the parent
% file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


theta = X(1);
beta = X(2);

E = T*Length(theta, beta, s) -F*Height(theta, beta, s);
