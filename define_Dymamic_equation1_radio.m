function[phi_z,Q_z,F_z]=define_Dymamic_equation1_radio(F_z,xvm_Nh,yvm_Nh,psi_h,P00_x,P00_y,sig_arw_0,sig_rrw_0,dt,k)
F_z(1,3) = -yvm_Nh(k-1);
F_z(2,3) = xvm_Nh(k-1);
phi_z = expm(F_z*dt);

q11(k-1) = ((cos(psi_h(k-1))).^2)*P00_x(1) + ((sin(psi_h(k-1))).^2)*P00_y(1);
q12(k-1) = sin(psi_h(k-1))*cos(psi_h(k-1))*P00_x(1)-sin(psi_h(k-1))*cos(psi_h(k-1))*P00_y(1);
q22(k-1) = ((sin(psi_h(k-1))).^2)*P00_x(1) + ((cos(psi_h(k-1))).^2)*P00_y(1);
q33(k-1) = sig_arw_0^2;
q44(k-1) = sig_rrw_0^2;

Q_z = [q11(k-1) q12(k-1)  0   0
       q12(k-1) q22(k-1)  0   0
        0   0  q33(k-1)  0
        0   0   0  q44(k-1)];

end