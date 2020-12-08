function [P00_z,K_z,z_update]=Kalman_Filter_update_RTK(P00_z,Pm_data,H,R,Pm_h,k)

zxm_z = Pm_data(k,1:2) - Pm_h(k,1:2); % 观测方程

Mu_z= zxm_z'; % u矩阵

z_update = zeros(8,1);

K_z=P00_z*H'/(H*P00_z*H'+R); % K gain
z_update = K_z*Mu_z; % z = K*u
%xz_h=xz_h+z_update;                        % x(k|k)
I=eye(8);
P00_z=(I-K_z*H)*P00_z;                     % P=(I-K*H)*P 是P矩阵更新公式

end

