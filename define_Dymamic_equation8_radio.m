function[phi_z,Q_z,F_z]=define_Dymamic_equation8_radio(F_z,Q_z,axm_h,aym_h,xvm_Nh,yvm_Nh,wzm_h,psi_h,sig_bx,sig_by,sig_xr,sig_yr,sig_arw_0,sig_rrw_0,dt,k,format_1)
%% format = 1 straight line ,format = 2 square  F为8*8矩阵公式参考PPT 这里选取 format_1==1 做计算
if (format_1 == 1)
F_z(2,3) = -cos(psi_h(k-1));
F_z(2,5) = -wzm_h(k-1);
F_z(2,6) = sin(psi_h(k-1));
F_z(2,7) = 1*(-sin(psi_h(k-1))*axm_h(k-1)-cos(psi_h(k-1))*aym_h(k-1));
F_z(2,8) = 1*(yvm_Nh(k-1));
F_z(5,2) = wzm_h(k-1);
F_z(5,3) = -sin(psi_h(k-1));
F_z(5,6) = -cos(psi_h(k-1));
F_z(5,7) = 1*(cos(psi_h(k-1))*axm_h(k-1)-sin(psi_h(k-1))*aym_h(k-1));
F_z(5,8) = 1*(-xvm_Nh(k-1));
elseif (format_1 == 2)
F_z(2,3) = 1*-cos(psi_h(k-1));
F_z(2,5) = -wzm_h(k-1);
F_z(2,6) = 1*sin(psi_h(k-1));
F_z(2,7) = 0*(-sin(psi_h(k-1))*axm_h(k-1)-cos(psi_h(k-1))*aym_h(k-1));
F_z(2,8) = 0*(yvm_Nh(k-1));
F_z(5,2) = wzm_h(k-1);
F_z(5,3) = 1*-sin(psi_h(k-1));
F_z(5,6) = 1*-cos(psi_h(k-1));
F_z(5,7) = 0*(cos(psi_h(k-1))*axm_h(k-1)-sin(psi_h(k-1))*aym_h(k-1));
F_z(5,8) = 0*(-xvm_Nh(k-1));
end
I=eye(8); % 8*8单位矩阵
%phi_z = expm(F_z*dt);
phi_z = I+(F_z*dt);  % phi为8*8矩阵
%phi_z = I + F_z*dt; 
Q_z(2,2) = ((cos(psi_h(k-1))).^2)*sig_bx^2 + ((sin(psi_h(k-1))).^2)*sig_by^2 + yvm_Nh(k-1)*yvm_Nh(k-1)*sig_arw_0^2; % Q为8*8矩阵参考PPT
Q_z(3,3) = sig_xr^2;
Q_z(5,5) = ((sin(psi_h(k-1))).^2)*sig_bx^2 + ((cos(psi_h(k-1))).^2)*sig_by^2 + xvm_Nh(k-1)*xvm_Nh(k-1)*sig_arw_0^2;
Q_z(6,6) = sig_yr^2;
Q_z(7,7) = 200*sig_arw_0^2;
Q_z(8,8) = 200*sig_rrw_0^2;

end