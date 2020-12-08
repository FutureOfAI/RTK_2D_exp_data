function [axm_h,aym_h,Vx_h,Vy_h,vxm_h,vym_h,bx_h,by_h]=initial_estimate_value1(m,vxm,vym,axm,aym,x_v,y_v,Vxerr,Vyerr)
%% 宣告估計參數
vxm_h=zeros(m);
vym_h=zeros(m);
axm_h=zeros(m);
aym_h=zeros(m);
Vx_h=zeros(m);
Vy_h=zeros(m);
by_h=zeros(m);
bx_h=zeros(m);

%% 設定初始值
vxm_h(1)=vxm(1);
vym_h(1)=vym(1);
axm_h(1)=axm(1);
aym_h(1)=aym(1);
Vx_h(1)=x_v(1) + Vxerr*rand;%x_v_hat(n) = x_v(n) + x1(n);
Vy_h(1)=y_v(1) + Vyerr*rand;%y_v_hat(n) = y_v(n) + x3(n);
by_h(1)=0;
bx_h(1)=0;
end