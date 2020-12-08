function [x_p,x_v,x_a,y_p,y_v,y_a] = trajectory1(radius,psi_0,wt,t0)

x_p = radius*sin(wt*t0)*cos(psi_0)+6.0;
x_v = radius*wt*cos(wt*t0)*cos(psi_0);
x_a = - radius*wt^2*sin(wt*t0)*cos(psi_0);
y_p = radius*sin(wt*t0)*sin(psi_0)+4.53;
y_v = radius*wt*cos(wt*t0)*sin(psi_0);
y_a = - radius*wt^2*sin(wt*t0)*sin(psi_0);
end
