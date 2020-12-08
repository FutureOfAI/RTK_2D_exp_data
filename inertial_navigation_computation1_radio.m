function[xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,wzm_h,psi_h]=inertial_navigation_computation1_radio(xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,Vx_h,Vy_h,wzm_h,psi_h,wzm,bz_h,k,dt)

wzm_h(k) = wzm(k) - bz_h(k-1);
wzm_h(k-1) = wzm(k-1) - bz_h(k-1);

psi_h(k) = psi_h(k-1) + (wzm_h(k-1)+wzm_h(k))*dt/2;
%psi_h(k) = psi_h(k-1) + wzm_h(k)*dt;
if (psi_h(k)>= 2*pi),
    psi_h(k) = psi_h(k) - 2*pi;
else
    psi_h(k) = psi_h(k);
end
%
xvm_Nh(k) = cos(psi_h(k))*Vx_h(k)-sin(psi_h(k))*Vy_h(k);
yvm_Nh(k) = sin(psi_h(k))*Vx_h(k)+cos(psi_h(k))*Vy_h(k);

xpm_Nh(k) = xpm_Nh(k-1)+ (xvm_Nh(k)+xvm_Nh(k-1))*dt/2;
ypm_Nh(k) = ypm_Nh(k-1)+ (yvm_Nh(k)+yvm_Nh(k-1))*dt/2;

end