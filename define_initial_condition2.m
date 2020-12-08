function [xz_h,P00_z]=define_initial_condition2(bz0,xperr,yperr,psierr,d2r)
%% 宣告狀態變數
%d2r=pi/180;
xz_h=zeros(4,1);      % x=[x1;x2;x3;x4;x5;x6;x7;x8]
%% 定義狀態之初始條件
xz_h(1,1) = xperr;            % x1(1),position error in x-axis (m)
xz_h(2,1) = yperr;          % x2(1),position error in y-axis (m)
xz_h(3,1) = psierr*d2r;            % x3(1),angle error in z-axis (rad)
xz_h(4,1) = bz0;          % x4(1),gyro bias error

P00_z = zeros(4,4);
P00_z(1,1) = xperr^2;
P00_z(2,2) = yperr^2;
P00_z(3,3) = (1*psierr*d2r)^2;
P00_z(4,4) =100* bz0^2;


end