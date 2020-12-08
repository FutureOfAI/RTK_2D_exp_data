function [P00_z,K_z,z_update]=Kalman_Filter_update4_radio(xz_h,P00_z,R1m,R2m,R3m,R4m,H,R,R1m_h,R2m_h,R3m_h,R4m_h,k)
%% ( k| k-1) ->( k| k)¡G­×¥¿ª¬ºA¥H¤ÎP
% x(:,k)=[x1(k);...;x8(k)],size of x = 8*n
%R1m_h(k) = sqrt((xr1-xpm_Nh(k))^2+(yr1-ypm_Nh(k))^2+h^2);
%R2m_h(k) = sqrt((xr2-xpm_Nh(k))^2+(yr2-ypm_Nh(k))^2+h^2);



zxm_z1(k) = R1m(k)-R1m_h(k)' ;
zxm_z2(k) = R2m(k)-R2m_h(k)' ;
zxm_z3(k) = R3m(k)-R3m_h(k)' ;
zxm_z4(k) = R4m(k)-R4m_h(k)' ;

%Mu_x=zxm_x(k)-zxm_hatx(k);
Mu_z=[zxm_z1(k)
    zxm_z2(k)
    zxm_z3(k)
    zxm_z4(k)];
K_z=P00_z*H'/(H*P00_z*H'+R);      % P(k|k-1)

%K(7) = 0;
%K(8) = 0;
z_update = K_z*Mu_z;
%xz_h=xz_h+z_update;                        % x(k|k)
I=eye(6);
P00_z=(I-K_z*H)*P00_z;                     % P(k|k)



end

