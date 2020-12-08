function[xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,psi_h,bz_h,bx_h,by_h]=upon_radiosensor_measurement_8_2(xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,k,psi_h,bz_h,bx_h,by_h,z_update)

xpm_Nh(k)=xpm_Nh(k)+z_update(1);
xvm_Nh(k)=xvm_Nh(k)+z_update(2);
bx_h(k) = bx_h(k) + z_update(3);
ypm_Nh(k)=ypm_Nh(k)+z_update(4);
yvm_Nh(k)=yvm_Nh(k)+z_update(5);
by_h(k) = by_h(k) + z_update(6);
psi_h(k)=psi_h(k)+z_update(7);%;
bz_h(k)=bz_h(k)+z_update(8);%;

end