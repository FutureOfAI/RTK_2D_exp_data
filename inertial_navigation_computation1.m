function [axm_h,aym_h,Vx_h,Vy_h]=inertial_navigation_computation1(axm_h,aym_h,Vx_h,Vy_h,axm,aym,bx_h,by_h,k,dt)
% k-1--> k :x_h,y_h,theta_h,Vx_h,Vy_h


axm_h(k)=axm(k)-bx_h(k);
aym_h(k)=aym(k)-by_h(k);

Vx_h(k)=Vx_h(k-1)+(axm_h(k-1)+axm_h(k))*dt/2;
Vy_h(k)=Vy_h(k-1)+(aym_h(k-1)+aym_h(k))*dt/2;
%
%Vx_h(k)=Vx_h(k-1)+axm_h(k-1)*dt;
%Vy_h(k)=Vy_h(k-1)+aym_h(k-1)*dt;


end