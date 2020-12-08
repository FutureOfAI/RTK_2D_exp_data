function [xz_h,P00_z]=define_initial_condition_8(bx0,by0,bz0,xperr,yperr,xverr,yverr,psierr,d2r)

xz_h=zeros(8,1);      

% xz_h(1,1) = 0*xperr;              % x1(1),position error in x-axis (m)
% xz_h(2,1) = 0*xverr;              % x2(1),velocity error in x-axis (m/sec)
xz_h(3,1) = 0*bx0;                % x3(1),accel bias error in x-axis (m/sec^2)
% xz_h(4,1) = 0*yperr;              % x4(1),position error in y-axis (m)
% xz_h(5,1) = 0*yverr;              % x5(1),velocity error in y-axis (m/sec)
xz_h(6,1) = 0*by0;                % x6(1),accel bias error in y-axis (m/sec^2)
% xz_h(7,1) = 0*psierr*d2r;         % x7(1),angle error in z-axis (rad)
xz_h(8,1) = bz0;       

P00_z = zeros(8,8);
P00_z(1,1) = xperr^2;
% P00_z(2,2) = 0*xverr^2;
P00_z(3,3) = 1*bx0^2;
P00_z(4,4) = yperr^2;
% P00_z(5,5) = 0*yverr^2;
P00_z(6,6) = 1*by0^2;
P00_z(7,7) = (1*psierr*d2r)^2;
P00_z(8,8) =100* bz0^2;


end