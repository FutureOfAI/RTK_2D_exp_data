function[xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,wzm_h,psi_h]=inertial_navigation_computation8_radio(xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,xam_Nh,yam_Nh,axm,aym,bx_h,by_h,psi_h,wzm,bz_h,k,dt,scale_factor_err)

axm_h(k)=axm(k)-bx_h(k);
aym_h(k)=aym(k)-by_h(k);
axm_h(k-1)=axm(k-1)-bx_h(k-1);
aym_h(k-1)=aym(k-1)-by_h(k-1);
% wzm_h(k-1) = (1-scale_factor_err)*wzm(k-1) - bz_h(k-1);
% wzm_h(k) = (1-scale_factor_err)*wzm(k) - bz_h(k-1);
wzm_h(k-1) = (1-scale_factor_err)*(wzm(k-1) - bz_h(k-1));
wzm_h(k) = (1-scale_factor_err)*(wzm(k) - bz_h(k-1));
% wzm_h(k-1) =wzm(k-1) - bz_h(k-1);
% wzm_h(k) = wzm(k) - bz_h(k-1);



psi_h(k) = psi_h(k-1) + (wzm_h(k-1)+wzm_h(k))*dt/2;
% psi_h(k) = psi_h(k-1) + wzm_h(k)*dt;
if (psi_h(k)>= 2*pi),
    psi_h(k) = psi_h(k) - 2*pi;
else
    psi_h(k) = psi_h(k);
end
%
% compute N-frame accelerations, velocities, and positions
xam_Nh(k) = cos(psi_h(k-1))*axm_h(k-1)-sin(psi_h(k-1))*aym_h(k-1)-wzm_h(k-1)*yvm_Nh(k-1);
yam_Nh(k) = sin(psi_h(k-1))*axm_h(k-1)+cos(psi_h(k-1))*aym_h(k-1)+wzm_h(k-1)*xvm_Nh(k-1);
xvm_Nh(k) = xvm_Nh(k-1)+(xam_Nh(k)+xam_Nh(k-1))*dt/2;
yvm_Nh(k) = yvm_Nh(k-1)+(yam_Nh(k)+yam_Nh(k-1))*dt/2;
xpm_Nh(k) = xpm_Nh(k-1)+ (xvm_Nh(k)+xvm_Nh(k-1))*dt/2;
ypm_Nh(k) = ypm_Nh(k-1)+ (yvm_Nh(k)+yvm_Nh(k-1))*dt/2;



end