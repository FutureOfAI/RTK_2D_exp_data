function[H,R,Pm_h]=RTK_discrete_EKF(xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,k,format_1)

Pm_h(k,1:2) = [xpm_Nh(k) ypm_Nh(k)];

r_partial_x = [1;0];
r_partial_y = [0;1];

H = [r_partial_x zeros(2,1) zeros(2,1) r_partial_y zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1)]; % HÎª8*2¾ØÕó

R = 1*[sig_x_r^2 0
    0 sig_y_r^2]; % R¾ØÕóÎª2*2

end