% ===============================================================
% Program name: flight_profile.m
% 7/20/2016
% ===============================================================
% Produce platform flight profile for calibrate the inertial sensor errors
% ===============================================================
clf;close all;clear;clc;
r2d=(180/pi);
d2r=(pi/180);
g = 9.8;                        % gravity constant in m/sec^2
% ===============================================================
% Generate roll angle and introduce bank angles to create yaw terms
% ===============================================================
% define the platform constant ground speed, roll angle acceleration
% ===============================================================
v = 10 + 0;                         % ground speed in m/sec
roll_a = 1.25*d2r;                 % band-band acceleration in m/sec^2
Ts = 4 + 0.0;
T1 = 14.3 + 0.0;
T_1 = Ts + T1;
T2 = 14.3;
T2 = T1;
T_2 = Ts + T2;
T3 = 14.3;
T3 = T2;
T_3 = Ts + T3;
T4 = 14.3;
T4 = T2;
T_4 = Ts + T4;
T5 = 14.3;
T5 = T2;
T_5 = Ts + T5;
T6 = 14.3;
T6 = T2;
T_6 = Ts + T6;
T7 = 14.3;
T7 = T2;
T_7 = Ts + T7;
T8 = 14.3;
T8 = T2;
T_8 = Ts + T8;
T = 4*Ts + 2*T1 + T2;
%
px0 = 54.958;
py0 = -189;
% ===============================================================
t = 0;
dt = 0.01;
N = (Ts+T1)/dt;
N2 = (Ts+T2)/dt;
%N = Ts/dt;
% ===============================================================
for i = 1:N,                % first 90 deg turn
    time(i) = t;
% generate roll angle
if (t <= Ts),
    if (t <= Ts/2),
        roll_angle(i) = (roll_a/2)*t^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(i) = -(roll_a/4)*Ts^2 + (roll_a*Ts)*t - (roll_a/2)*t^2;
    end
else
    roll_angle(i) = roll_angle(i-1);
end
% generate yaw angle rate
yaw_rate(i) = (g*tan(roll_angle(i)))/v;
% produce the yaw angle by integrating the above yaw rate
if (i == 1),
    yaw_angle(i) = 0;
else
    yaw_angle(i) = yaw_angle(i-1) + yaw_rate(i)*dt;
end
% compute platform vecolities
vx(i) = v*cos(yaw_angle(i));
vy(i) = v*sin(yaw_angle(i));
ax(i) = -v*yaw_rate(i)*sin(yaw_angle(i));
ay(i) = v*yaw_rate(i)*cos(yaw_angle(i));
%
f27(i) = -ay(i) + yaw_rate(i)*vx(i);
f57(i) = ax(i) + yaw_rate(i)*vy(i);
% compute platform positions
if (i ==1),
    px(i) = px0;
    py(i) = py0;
else
px(i) = px(i-1) + vx(i)*dt;
py(i) = py(i-1) + vy(i)*dt;
end
    t = t + dt;
end
% ============================================
% end of first roll or bank turn (positive roll turn)
% ============================================
% start the next bank turn (negative bank turn)
% ============================================
%roll_a = - roll_a;
for k = i+1:i+N,
   time(k) = t;
   tt = t - T_1;
% generate roll angle
if (t <= T_1+Ts),
    if (t <= T_1+Ts/2),
        roll_angle(k) = roll_angle(i) - (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k) = (roll_a/2)*Ts^2 - (roll_a*Ts)*tt + (roll_a/2)*tt^2;
    end
else
    roll_angle(k) = roll_angle(k-1);
end
% generate yaw angle rate
yaw_rate(k) = (g*tan(roll_angle(k)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k == 1),
    yaw_angle(k) = 0;
else
    yaw_angle(k) = yaw_angle(k-1) + yaw_rate(k)*dt;
end
% compute platform vecolities
vx(k) = v*cos(yaw_angle(k));
vy(k) = v*sin(yaw_angle(k));
ax(k) = -v*yaw_rate(k)*sin(yaw_angle(k));
ay(k) = v*yaw_rate(k)*cos(yaw_angle(k));
%
f27(k) = -ay(k) + yaw_rate(k)*vx(k);
f57(k) = ax(k) + yaw_rate(k)*vy(k);
% compute platform positions
if (k ==1),
    px(k) = px0;
    py(k) = py0;
else
px(k) = px(k-1) + vx(k)*dt;
py(k) = py(k-1) + vy(k)*dt;
end
    t = t + dt;
end
% =========================================================
% end of second roll or bank turn 
% =========================================================
% Ready for next bank turn (positive band turn again)
% =========================================================
for k1 = k+1:k+N,
   time(k1) = t;
   tt = t - T_1 - T_2;
% generate roll angle
if (t <= T_1+T_2+Ts),
    if (t <= T_1+T_2+Ts/2),
        roll_angle(k1) = roll_angle(k) + (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k1) = -(roll_a/4)*Ts^2 + (roll_a*Ts)*tt - (roll_a/2)*tt^2;
    end
else
    roll_angle(k1) = roll_angle(k1-1);
end
% generate yaw angle rate
yaw_rate(k1) = (g*tan(roll_angle(k1)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k1 == 1),
    yaw_angle(k1) = 0;
else
    yaw_angle(k1) = yaw_angle(k1-1) + yaw_rate(k1)*dt;
end
% compute platform vecolities
vx(k1) = v*cos(yaw_angle(k1));
vy(k1) = v*sin(yaw_angle(k1));
ax(k1) = -v*yaw_rate(k1)*sin(yaw_angle(k1));
ay(k1) = v*yaw_rate(k1)*cos(yaw_angle(k1));
%
f27(k1) = -ay(k1) + yaw_rate(k1)*vx(k1);
f57(k1) = ax(k1) + yaw_rate(k1)*vy(k1);
% compute platform positions
if (k1 ==1),
    px(k1) = px0;
    py(k1) = py0;
else
px(k1) = px(k1-1) + vx(k1)*dt;
py(k1) = py(k1-1) + vy(k1)*dt;
end
    t = t + dt;
end
% =========================================================
% end of third roll or bank turn 
% =========================================================
% Ready for next bank turn (negative band turn again)
% =========================================================
for k2 = k1+1:k1+N,
   time(k2) = t;
   tt = t - T_1 - T_2 - T_3;
% generate roll angle
if (t <= T_1+T_2+T_3+Ts),
    if (t <= T_1+T_2+T_3+Ts/2),
        roll_angle(k2) = roll_angle(k1) - (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k2) = (roll_a/2)*Ts^2 - (roll_a*Ts)*tt + (roll_a/2)*tt^2;
    end
else
    roll_angle(k2) = roll_angle(k2-1);
end
% generate yaw angle rate
yaw_rate(k2) = (g*tan(roll_angle(k2)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k2 == 1),
    yaw_angle(k2) = 0;
else
    yaw_angle(k2) = yaw_angle(k2-1) + yaw_rate(k2)*dt;
end
% compute platform vecolities
vx(k2) = v*cos(yaw_angle(k2));
vy(k2) = v*sin(yaw_angle(k2));
ax(k2) = -v*yaw_rate(k2)*sin(yaw_angle(k2));
ay(k2) = v*yaw_rate(k2)*cos(yaw_angle(k2));
%
f27(k2) = -ay(k2) + yaw_rate(k2)*vx(k2);
f57(k2) = ax(k2) + yaw_rate(k2)*vy(k2);
% compute platform positions
if (k2 ==1),
    px(k2) = px0;
    py(k2) = py0;
else
px(k2) = px(k2-1) + vx(k2)*dt;
py(k2) = py(k2-1) + vy(k2)*dt;
end
    t = t + dt;
end
% =========================================
% end of fouth roll or bank turn 
% =========================================
% Ready for next bank turn (positive band turn again)
% ==================================================
for k3 = k2+1:k2+N,
   time(k3) = t;
   tt = t - T_1 - T_2 - T_3 - T_4;
% generate roll angle
if (t <= T_1+T_2+T_3+T_4+Ts),
    if (t <= T_1+T_2+T_3++T_4+Ts/2),
        roll_angle(k3) = roll_angle(k2) + (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k3) = -(roll_a/4)*Ts^2 + (roll_a*Ts)*tt - (roll_a/2)*tt^2;
    end
else
    roll_angle(k3) = roll_angle(k3-1);
end
% generate yaw angle rate
yaw_rate(k3) = (g*tan(roll_angle(k3)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k3 == 1),
    yaw_angle(k3) = 0;
else
    yaw_angle(k3) = yaw_angle(k3-1) + yaw_rate(k3)*dt;
end
% compute platform vecolities
vx(k3) = v*cos(yaw_angle(k3));
vy(k3) = v*sin(yaw_angle(k3));
ax(k3) = -v*yaw_rate(k3)*sin(yaw_angle(k3));
ay(k3) = v*yaw_rate(k3)*cos(yaw_angle(k3));
%
f27(k3) = -ay(k3) + yaw_rate(k3)*vx(k3);
f57(k3) = ax(k3) + yaw_rate(k3)*vy(k3);
% compute platform positions
if (k3 ==1),
    px(k3) = px0;
    py(k3) = py0;
else
px(k3) = px(k3-1) + vx(k3)*dt;
py(k3) = py(k3-1) + vy(k3)*dt;
end
    t = t + dt;
end
% =========================================
% end of fifth roll or bank turn 
% =========================================
% Ready for next bank turn (negative band turn again)
% ===============================================
for k4 = k3+1:k3+N,
   time(k4) = t;
   tt = t - T_1 - T_2 - T_3 - T_4 - T_5;
% generate roll angle
if (t <= T_1+T_2+T_3+T_4+T_5+Ts),
    if (t <= T_1+T_2+T_3+T_4+T_5+Ts/2),
        roll_angle(k4) = roll_angle(k3) - (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k4) = (roll_a/2)*Ts^2 - (roll_a*Ts)*tt + (roll_a/2)*tt^2;
    end
else
    roll_angle(k4) = roll_angle(k4-1);
end
% generate yaw angle rate
yaw_rate(k4) = (g*tan(roll_angle(k4)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k4 == 1),
    yaw_angle(k4) = 0;
else
    yaw_angle(k4) = yaw_angle(k4-1) + yaw_rate(k4)*dt;
end
% compute platform vecolities
vx(k4) = v*cos(yaw_angle(k4));
vy(k4) = v*sin(yaw_angle(k4));
ax(k4) = -v*yaw_rate(k4)*sin(yaw_angle(k4));
ay(k4) = v*yaw_rate(k4)*cos(yaw_angle(k4));
%
f27(k4) = -ay(k4) + yaw_rate(k4)*vx(k4);
f57(k4) = ax(k4) + yaw_rate(k4)*vy(k4);
% compute platform positions
if (k4 ==1),
    px(k4) = px0;
    py(k4) = py0;
else
px(k4) = px(k4-1) + vx(k4)*dt;
py(k4) = py(k4-1) + vy(k4)*dt;
end
    t = t + dt;
end
% =========================================
% end of six roll or bank turn 
% ========================================
% Ready for next bank turn (positive band turn again)
% ==================================================
for k5 = k4+1:k4+N,
   time(k5) = t;
   tt = t - T_1 - T_2 - T_3 - T_4 - T_5 - T_6;
% generate roll angle
if (t <= T_1+T_2+T_3+T_4+T_5+T_6+Ts),
    if (t <= T_1+T_2+T_3++T_4+T_5+T_6+Ts/2),
        roll_angle(k5) = roll_angle(k4) + (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k5) = -(roll_a/4)*Ts^2 + (roll_a*Ts)*tt - (roll_a/2)*tt^2;
    end
else
    roll_angle(k5) = roll_angle(k5-1);
end
% generate yaw angle rate
yaw_rate(k5) = (g*tan(roll_angle(k5)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k5 == 1),
    yaw_angle(k5) = 0;
else
    yaw_angle(k5) = yaw_angle(k5-1) + yaw_rate(k5)*dt;
end
% compute platform vecolities
vx(k5) = v*cos(yaw_angle(k5));
vy(k5) = v*sin(yaw_angle(k5));
ax(k5) = -v*yaw_rate(k5)*sin(yaw_angle(k5));
ay(k5) = v*yaw_rate(k5)*cos(yaw_angle(k5));
%
f27(k5) = -ay(k5) + yaw_rate(k5)*vx(k5);
f57(k5) = ax(k5) + yaw_rate(k5)*vy(k5);
% compute platform positions
if (k5 ==1),
    px(k5) = px0;
    py(k5) = py0;
else
px(k5) = px(k5-1) + vx(k5)*dt;
py(k5) = py(k5-1) + vy(k5)*dt;
end
    t = t + dt;
end
% =========================================
% end of seventh roll or bank turn 
% =========================================
% Ready for next bank turn (negative band turn again)
% ===============================================
for k6 = k5+1:k5+N,
   time(k6) = t;
   tt = t - T_1 - T_2 - T_3 - T_4 - T_5 - T_6 - T_7;
% generate roll angle
if (t <= T_1+T_2+T_3+T_4+T_5+T_6+T_7+Ts),
    if (t <= T_1+T_2+T_3+T_4+T_5+T_6+T_7+Ts/2),
        roll_angle(k6) = roll_angle(k5) - (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k6) = (roll_a/2)*Ts^2 - (roll_a*Ts)*tt + (roll_a/2)*tt^2;
    end
else
    roll_angle(k6) = roll_angle(k6-1);
end
% generate yaw angle rate
yaw_rate(k6) = (g*tan(roll_angle(k6)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k6 == 1),
    yaw_angle(k6) = 0;
else
    yaw_angle(k6) = yaw_angle(k6-1) + yaw_rate(k6)*dt;
end
% compute platform vecolities
vx(k6) = v*cos(yaw_angle(k6));
vy(k6) = v*sin(yaw_angle(k6));
ax(k6) = -v*yaw_rate(k6)*sin(yaw_angle(k6));
ay(k6) = v*yaw_rate(k6)*cos(yaw_angle(k6));
%
f27(k6) = -ay(k6) + yaw_rate(k6)*vx(k6);
f57(k6) = ax(k6) + yaw_rate(k6)*vy(k6);
% compute platform positions
if (k6 ==1),
    px(k6) = px0;
    py(k6) = py0;
else
px(k6) = px(k6-1) + vx(k6)*dt;
py(k6) = py(k6-1) + vy(k6)*dt;
end
    t = t + dt;
end
% =========================================
% end of eighth roll or bank turn 
% =========================================
% Ready for next bank turn (positive band turn again)
% ==================================================
for k7 = k6+1:k6+N,
   time(k7) = t;
   tt = t - T_1 - T_2 - T_3 - T_4 - T_5 - T_6 - T_7 - T_8;
% generate roll angle
if (t <= T_1+T_2+T_3+T_4+T_5+T_6+T_7+T_8+Ts),
    if (t <= T_1+T_2+T_3++T_4+T_5+T_6+T_7+T_8+Ts/2),
        roll_angle(k7) = roll_angle(k6) + (roll_a/2)*tt^2;
    else
        %roll_angle(i) = (roll_a/2)*Ts^2 - (5*roll_a*Ts/4)*t + (roll_a)*t^2;
        roll_angle(k7) = -(roll_a/4)*Ts^2 + (roll_a*Ts)*tt - (roll_a/2)*tt^2;
    end
else
    roll_angle(k7) = roll_angle(k7-1);
end
% generate yaw angle rate
yaw_rate(k7) = (g*tan(roll_angle(k7)))/v;
% produce the yaw angle by integrating the above yaw rate
if (k5 == 1),
    yaw_angle(k7) = 0;
else
    yaw_angle(k7) = yaw_angle(k7-1) + yaw_rate(k7)*dt;
end
% compute platform vecolities
vx(k7) = v*cos(yaw_angle(k7));
vy(k7) = v*sin(yaw_angle(k7));
ax(k7) = -v*yaw_rate(k7)*sin(yaw_angle(k7));
ay(k7) = v*yaw_rate(k7)*cos(yaw_angle(k7));
%
f27(k7) = -ay(k7) + yaw_rate(k7)*vx(k7);
f57(k7) = ax(k7) + yaw_rate(k7)*vy(k7);
% compute platform positions
if (k7 ==1),
    px(k7) = px0;
    py(k7) = py0;
else
px(k7) = px(k7-1) + vx(k7)*dt;
py(k7) = py(k7-1) + vy(k7)*dt;
end
    t = t + dt;
end
% =========================================
% end of nineth roll or bank turn 
% =========================================
%
figure (1)
plot(time,roll_angle*r2d)
xlabel('Time in seconds')
ylabel('Roll angle in deg')
grid
%
figure (2)
subplot (211)
plot(time,yaw_rate*r2d)
%xlabel('Time in seconds')
ylabel('yaw angle rate in deg/sec')
grid
subplot(212)
plot(time,vx,'r',time,vy,'b')
ylabel('Velocity in m/sec')
grid
%
figure (3)
plot(time,yaw_angle*r2d)
xlabel('Time in seconds')
ylabel('yaw angle in deg')
grid
%
figure (4)
subplot(211)
plot(time,f27)
ylabel('Coef f27 in m/sec^2')
grid
subplot(212)
plot(time,f57)
xlabel('Time in seconds')
ylabel('Coef f57 in m/sec^2')
grid
%
figure (5)
plot(px,py)
xlabel('Platform x-axis position in m')
ylabel('Platform y-axis position in m')
grid
%
figure (6)
subplot(211)
plot(time,roll_angle*r2d)
ylabel('Platform roll angle in deg')
grid
subplot(212)
plot(time,yaw_angle*r2d)
xlabel('Time in seconds')
ylabel('Platform yaw angle in deg')
grid
%
figure (7)
subplot(211)
plot(time,ax/g)
ylabel('X accel in g')
grid
subplot(212)
plot(time,ay/g)
xlabel('Time in seconds')
ylabel('Y accel in g')
grid

        
