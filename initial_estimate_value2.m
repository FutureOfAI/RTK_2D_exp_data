function [xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,wzm_h,psi_h,bz_h]=initial_estimate_value2(m,x_v_N,y_v_N,x_p_N,y_p_N,wzm,psi,d2r,xverr,yverr,xperr,yperr,psierr)
%4-state Kalman filter
xvm_Nh=zeros(m);
yvm_Nh=zeros(m);
xpm_Nh=zeros(m);
ypm_Nh=zeros(m);
wzm_h=zeros(m);
psi_h=zeros(m);
bz_h=zeros(m);

xvm_Nh(1)=x_v_N(1) + xverr*rand;
yvm_Nh(1)=y_v_N(1) + yverr*rand;
xpm_Nh(1)=x_p_N(1) + xperr*rand;
ypm_Nh(1)=y_p_N(1) + yperr*rand;
wzm_h(1)=wzm(1);
psi_h(1)=psi(1) + psierr*rand*d2r;
bz_h(1)=0;
end