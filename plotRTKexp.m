close all;

n1 = 400;
%n2 = 100;
n2 = k-1;

figure (1)
plot(Pm_IMU_h(1:n2,1),Pm_IMU_h(1:n2,2),'b.',xpm_Nh(1:n2),ypm_Nh(1:n2),'r--','LineWidth',2)
hold on
plot(Pm_IMU_h(2,1),Pm_IMU_h(2,2),'diamond',Pm_IMU_h(n2,1),Pm_IMU_h(n2,2),'rsquare','MarkerSize',13,'LineWidth',2)%,5.4,0.74,'bx'
legend('RTK path','EKF path','Start point','End point');
xlabel('X position in m')
ylabel('Y position in m')
axis([2 9 -2 7]);
set(gca,'YDir','reverse')
title('EKF RTK random positioning')
grid
%
figure (2)
subplot(311)
plot(t00,accx(1:k)/g,'r')
xlabel('Time in sec')
ylabel('X accel in g')
subplot(312)
plot(t00,accy(1:k)/g,'r')
xlabel('Time in sec')
ylabel('Y accel in g')
subplot(313)
plot(t00,groz(1:k),'r')
xlabel('Time in sec')
ylabel('Z gyro in rad/s')
%
figure(3);
subplot(3,1,1), plot(t0,bx_h/g,'r');         
xlabel('Time in seconds');ylabel('Accel bias along-axis in g');
subplot(3,1,2), plot(t0,by_h/g,'r');         
xlabel('Time in seconds');ylabel('Accel bias perp-axis in g');
subplot(3,1,3), plot(t0,bz_h*r2d,'r');
xlabel('Time in seconds');ylabel('Gyro bias in deg');
%
% figure (4)
% subplot(3,1,1), plot(t0,bx/g,'r');         
% xlabel('Time in seconds');ylabel('Accel bias along-axis in g');
% subplot(3,1,2), plot(t0,by/g,'r');         
% xlabel('Time in seconds');ylabel('Accel bias perp-axis in g');
% subplot(3,1,3), plot(t0,bz*r2d,'r');
% xlabel('Time in seconds');ylabel('Gyro bias in deg/sec');
%
% figure (5)
% subplot(2,1,1), plot(t0(n1:n2),(bx(n1:n2)-bx_h(n1:n2))/g,'r');         
% xlabel('Time in seconds');ylabel('Accel bias err in x-axis in g');
% grid
% subplot(2,1,2), plot(t0(n1:n2),(by(n1:n2)-by_h(n1:n2))/g,'r');         
% xlabel('Time in seconds');ylabel('Accel bias err in y-axis in g');
% grid
%
figure (6)
plot(t0,(psi_h')*r2d)
xlabel('Time in seconds');ylabel('psi_h Angle in deg');
%
figure (7)
subplot(211)
plot(100:k,Pm_IMU_h(100:k,1)-xpm_Nh(100:k),'r')
xlabel('Time in sec')
ylabel('x postion err in m')
subplot(212)
plot(100:k,Pm_IMU_h(100:k,2)-ypm_Nh(100:k),'r')
xlabel('Time in sec')
ylabel('y postion err in m')
%
figure (8)
% subplot(211)
plot(100:k,psi_h(100:k)*r2d,'r.',100:k,attz(100:k)*r2d,'g.')
legend('EKF fusion angle','IMU module angle');
xlabel('Time in sec')
ylabel('Euler yaw in deg')
% subplot(212)
% plot(100:k,(psi_h(100:k)-attz(100:k))*r2d,'r')
% xlabel('Time in sec')
% ylabel('yaw err in deg')
%
figure (9)
subplot(211)
plot(100:k,xvm_Nh(100:k),'r')
xlabel('Time in sec')
ylabel('x velocity in m/s')
subplot(212)
plot(100:k,yvm_Nh(100:k),'r')
xlabel('Time in sec')
ylabel('y velocity in m/s')
%

