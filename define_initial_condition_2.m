function [xz_h,P00_z]=define_initial_condition_2(bz0,xperr,yperr,xverr,yverr,psierr,d2r)

xz_h=zeros(6,1);      
%% 定義狀態之初始條件
xz_h(1,1) = xperr;              % x1(1),position error in x-axis (m)
xz_h(2,1) = xverr;              % x2(1),velocity error in x-axis (m/sec)
xz_h(3,1) = yperr;              % x3(1),position error in y-axis (m)
xz_h(4,1) = yverr;              % x4(1),velocity error in y-axis (m/sec)
xz_h(5,1) = psierr*d2r;         % x5(1),angle error in z-axis (rad)
xz_h(6,1) = bz0;                % x6(1),gyro bias error

P00_z = zeros(6,6);
P00_z(1,1) = xperr^2;
P00_z(2,2) = xverr^2;
P00_z(3,3) = yperr^2;
P00_z(4,4) = yverr^2;
P00_z(5,5) = (1*psierr*d2r)^2;
P00_z(6,6) =100* bz0^2;


end