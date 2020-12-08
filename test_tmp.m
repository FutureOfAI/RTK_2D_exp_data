close all;
clear all;
clc;
% load experiment data
gro_Data = csvread('10_State_data0202.csv');
r2d = 180/pi;
dt = 0.01;
angle = zeros(2,100);
for i=1:2
    for j = 2:100
        angle(i,j) = angle(i,j-1) + gro_Data(i,j)*dt;
    end
end

mean(angle(1,93:98)*r2d)
mean(angle(2,90:94)*r2d)


timer = 1:100;
figure (1)
subplot(211)
plot(timer,gro_Data(1,timer))
ylabel('X gyro in rad')
grid
subplot(212)
plot(timer,gro_Data(2,timer))
ylabel('Y gyro in rad')
grid
%
figure (2)
subplot(211)
plot(timer,angle(1,timer)*r2d)
ylabel('po angle in rad')
grid
subplot(212)
plot(timer,angle(2,timer)*r2d)
ylabel('ne angle in rad')
grid