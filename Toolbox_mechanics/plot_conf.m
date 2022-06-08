function plot_conf(theta0, phi0, s, str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the sphape of the vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = s/(sin(theta0)+cos(phi0));

ang = linspace(-pi/2+theta0,pi-phi0);
plot(R*cos(ang)-min(R*cos(ang)), R*sin(ang)-min(R*sin(ang)), str,'LineWidth',3)
hold on
plot(-R*cos(ang)+min(R*cos(ang)), R*sin(ang)-min(R*sin(ang)), str,'LineWidth',3)
axis equal
