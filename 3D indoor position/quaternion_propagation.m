% quaternion_propagation.m
% scripts to perform quaternion propagation using third-order approximation
% other subrountines needed
% cx, cy, cz, qmult, q2dc, DTQ
% 4/24/2017
% Andy Wu
% ============================================================================
clear
format short e
clc
% ==========================================================================
% define some constants
rtd = 180/pi;                       % radian to degree conversion factor                                       
dtr = 1/rtd;                        % degree to radian
atr = (1/3600)*dtr;                 % arc sec to radian
orbit_rate = 2*pi/(1*10);           % bus orbit rate in rad/sec - 1 min LEO robit
period = 2*pi/orbit_rate;
t = 0;
dt = 1/100;                         % time period for attitude propagation
N = round(period/dt);
N = 2*N;
% =========================================================================
% define three initial Euler angle
% =========================================================================
roll_0 = 30*dtr;
pitch_0 = 0*dtr;
yaw_0 = 0*dtr;
roll_rate = 0.1*dtr;
pitch_rate = 0.1*dtr;
yaw_rate = 2*pi*(1/10);
r1 = [1 0 0];
r2 = [0 1 0];
r3 = [0 0 1];
% =========================================================================
% Generate w_EB_B
% =========================================================================
for n1 = 1:N,
    time(n1) = t;
    roll(n1) = roll_0 + roll_rate*t;
    pitch(n1) = pitch_0 + pitch_rate*t;
    yaw(n1) = yaw_0 + yaw_rate*t;
% generate inertial rates
% =========================================================
    T_1 = r1'*r1 + inv(cx(roll(n1)))*(r2'*r2) + inv(cx(roll(n1)))*inv(cy(pitch(n1)))*(r3'*r3);
    w_EB_B(1:3,n1) = T_1*[roll_rate pitch_rate yaw_rate]';
% =========================================================
    t = t + dt;
end
% =========================================================================
figure (1)
plot (time,w_EB_B)
xlabel('Time in seconds')
ylabel('Rates in rad.')
grid
% =========================================================================
% Quaternion propagation
% =========================================================================
CE_B = cz(yaw(1))*cy(pitch(1))*cx(roll(1));
QE_B = DTQ(CE_B);
QE_B_m = QE_B;
QB_E_m = -QE_B_m;
QB_E_m(4,1) = QE_B_m(4,1);
dQ = [0 0 0 1]';
dq1(1) = 2*dQ(1,1);
dq2(2) = 2*dQ(2,1);
dq3(1) = 2*dQ(3,1);
for n = 2:N,
    CE_B = cz(yaw(n))*cy(pitch(n))*cx(roll(n));
    QE_B = DTQ(CE_B);
    QB_E = -QE_B;
    QB_E(4,1) = QE_B(4,1);
% ============================================
% quaternion propagation algorithm
% ============================================
    d1 = w_EB_B(1,n)*dt/2;
    d2 = w_EB_B(2,n)*dt/2;
    d3 = w_EB_B(3,n)*dt/2;
    d1p = w_EB_B(1,n-1)*dt/2;
    d2p = w_EB_B(2,n-1)*dt/2;
    d3p = w_EB_B(3,n-1)*dt/2;
    d0_s = d1^2 + d2^2 + d3^2;
    q1 = d1 - (d0_s*d1 + d3p*d2 - d2p*d3)/6;
    q2 = d2 - (d0_s*d2 + d1p*d3 - d3p*d1)/6;
    q3 = d3 - (d0_s*d3 + d2p*d1 - d1p*d2)/6;
    q4 = 1 - d0_s/2;
    delta_Q = [q1 q2 q3 q4]';
    QE_B_m = qmult(delta_Q,QE_B_m);
% =============================================
% check the angle errors
% =============================================
    inv_QE_B_m = -QE_B_m;
    inv_QE_B_m(4,1) = QE_B_m(4,1);
    dQ1 = qmult(QE_B,inv_QE_B_m);
    dq11(n) = 2*dQ1(1,1);
    dq21(n) = 2*dQ1(2,1);
    dq31(n) = 2*dQ1(3,1);
end
% ============================================
figure (2)
plot (time,dq11*1e06,'b',time,dq21*1e06,'g',time,dq31*1e06,'r')
xlabel('Time in seconds')
ylabel('Angle errors in microrad')
grid

