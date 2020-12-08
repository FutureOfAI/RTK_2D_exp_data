function [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,by_h,bx_h,wzm_h,psi_h,bz_h]=initial_estimate_value8_radio(m,x_a_N,y_a_N,x_v_N,y_v_N,x_p_N,y_p_N,wzm,axm,aym,psi,d2r,xverr,yverr,xperr,yperr,xaerr,yaerr,psierr,position_x,position_y,bz0)
%
axm_h=zeros(m);
aym_h=zeros(m);
by_h=zeros(m);
bx_h=zeros(m);

xvm_Nh=zeros(m);
yvm_Nh=zeros(m);
xpm_Nh=zeros(m);
ypm_Nh=zeros(m);
xam_Nh=zeros(m);
yam_Nh=zeros(m);
wzm_h=zeros(m);
psi_h=zeros(m);
bz_h=zeros(m);

axm_h(1)=axm(1);
aym_h(1)=aym(1);
by_h(1)=0*9.8;
bx_h(1)=0*9.8;

xvm_Nh(1)=0;%14s-1.5609;%2.2313%0*x_v_N(1) ;%+ xverr*rand;
yvm_Nh(1)=0;%14s-1.4460;%1.5576%0*y_v_N(1); %+ yverr*rand;
% xpm_Nh(1)=4.7058;%4.8;%14s4.538;%x_p_N(1); %+ xperr*rand;
% ypm_Nh(1)=6.7240;%6.4;%14s6.066;%y_p_N(1) ;% + yperr*rand;  aj_1
xpm_Nh(1)=position_x;%4.7158;%4.8;%14s4.538;%x_p_N(1); %+ xperr*rand;
ypm_Nh(1)=position_y;%6.6467;%6.4;%14s6.066;%y_p_N(1) ;% + yperr*rand;
xam_Nh(1)=0*x_a_N(1); %+ xaerr*rand;
yam_Nh(1)=0*y_a_N(1); %+ yaerr*rand;
wzm_h(1)=wzm(1);
psi_h(1)=-5*d2r;%psi(1) ;%+ psierr*rand*d2r;-180
bz_h(1)=bz0;
% bz_h(1)=-4.37*d2r;

end