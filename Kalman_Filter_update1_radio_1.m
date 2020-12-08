function [xz_h,P00_z,K_z,z_update]=Kalman_Filter_update1_radio_1(xz_h,P00_z,x_p_N,y_p_N,xm,ym,xpm_Nh,ypm_Nh,H,R,k)
%% ( k| k-1) ->( k| k)¡G­×¥¿ª¬ºA¥H¤ÎP
% x(:,k)=[x1(k);...;x8(k)],size of x = 8*n
%R1m_h(k) = sqrt((xr1-xpm_Nh(k))^2+(yr1-ypm_Nh(k))^2+h^2);
%R2m_h(k) = sqrt((xr2-xpm_Nh(k))^2+(yr2-ypm_Nh(k))^2+h^2);
x_p_Nm = x_p_N(k) + xm*randn;
y_p_Nm = y_p_N(k) + ym*randn;

zxm_z1(k) = x_p_Nm - xpm_Nh(k) ;
zxm_z2(k) = y_p_Nm - ypm_Nh(k) ;

%Mu_x=zxm_x(k)-zxm_hatx(k);
Mu_z=[zxm_z1(k);zxm_z2(k)];
K_z=P00_z*H'/(H*P00_z*H'+R);      % P(k|k-1)

%K(7) = 0;
%K(8) = 0;
z_update = K_z*Mu_z;
xz_h=xz_h+z_update;                        % x(k|k)
I=eye(4);
P00_z=(I-K_z*H)*P00_z;                     % P(k|k)



end

