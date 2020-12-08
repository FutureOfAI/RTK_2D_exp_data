function[phi_z,Q_z,F_z]=define_Dymamic_equation2_radio(F_z,Q_z,axm_h,aym_h,xvm_Nh,yvm_Nh,wzm_h,psi_h,P00_x,P00_y,sig_arw_0,sig_rrw_0,dt,k)
F_z(2,4) = -wzm_h(k-1);
F_z(2,5) = -(sin(psi_h(k-1))*axm_h(k-1)+cos(psi_h(k-1))*aym_h(k-1));
F_z(2,6) = yvm_Nh(k-1);
F_z(4,2) = wzm_h(k-1);
F_z(4,5) = (cos(psi_h(k-1))*axm_h(k-1)-sin(psi_h(k-1))*aym_h(k-1));
F_z(4,6) = -xvm_Nh(k-1);
phi_z = expm(F_z*dt);

Q_z(5,5) = sig_arw_0^2;
Q_z(6,6) = sig_rrw_0^2;
Q_z(2,2) = ((cos(psi_h(k-1))).^2)*P00_x(4) + ((sin(psi_h(k-1))).^2)*P00_y(4)+(yvm_Nh(k-1))^2*Q_z(5,5);
Q_z(4,4) = ((sin(psi_h(k-1))).^2)*P00_x(4) + ((cos(psi_h(k-1))).^2)*P00_y(4)+(xvm_Nh(k-1))^2*Q_z(5,5);


end