function [P00_z,K_z,z_update]=Kalman_Filter_update_mult_radio(P00_z,Rm_data,H,R,Rm_h,k)

zxm_z = Rm_data(k,1:4) - Rm_h(k,1:4);

Mu_z= zxm_z';

z_update = zeros(8,1);

K_z=P00_z*H'/(H*P00_z*H'+R);    
z_update = K_z*Mu_z;
%xz_h=xz_h+z_update;                        % x(k|k)
I=eye(8);
P00_z=(I-K_z*H)*P00_z;                     % P(k|k)

end

