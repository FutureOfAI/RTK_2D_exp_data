function[x_h,y_h,vxm_h_N,vym_h_N,xm_h_N,ym_h_N,wzm_h,psi_h,bx_h,by_h=initial_estimate_value_radio(m,vxm,vym,Vx_h,Vy_h,wzm)
%% 宣告估計參數
x_h=zeros(m);
y_h=zeros(m);
vxm_h_N=zeros(m);
vym_h_N=zeros(m);
xm_h_N=zeros(m);
ym_h_N=zeros(m);
wzm_h=zeros(m);
psi_h=zeros(m);
bx_h=zeros(m);
by_h=zeros(m);

%% 設定初始值
x_h(1)=xm(1);
y_h(1)=ym(1);
theta_h(1)=thetam(1);
wzm_h(1)=wzm(1);
aam_h(1)=aam(1);
apm_h(1)=apm(1);
Vx_h(1)=0 + 0.*x_dot(1);
Vy_h(1)=0 + 0.*y_dot(1);
bg_h(1)=0;
bp_h(1)=0;
ba_h(1)=0;
end