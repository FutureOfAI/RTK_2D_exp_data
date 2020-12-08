% Program name: radio_sensor_8_2.m
% A new program that deal with non-constant psi angle where
% a 8-state Kalman filter with only 2 radio sensors without optical flow
% sensors
close all;
clear all;
clc;
format short e
r2d = (180/pi);
d2r = (pi/180);
g =9.8;
scale_factor_err = -0.0512;%; °f -0.0512 ¶¶ -0.2577-0.06112
format_1 = 1; % stiaight line = 1 (-0.89426), square = 2 (-4.37)  
%
for ii = 1:100,
%======================================================
% Design parameters to be tuned
%=====================================================
% Select motion profiles
% trajectory is generated in navigation frame (N-frame)
%================================
% Set motion profile flags
% profile_flag = 1:         straight motion
% profile_flag = 2;         circular motion
% profile_flag = 3;         race track motion
profile_flag = 1;
%
if (profile_flag ==1),
% =====================================================
% Straight motion with psi = constant angle, or wn = 0
% =====================================================
% positions of two ratio sensors's positions
% ================================
xr1 = 0;% in meter
yr1 = 25;% in meter
xr2 = 0;
yr2 = -25;
xr3 = 25;% in meter
yr3 = 0;% in meter
xr4 = -25;
yr4 = 0;
h = 4;% ceiling high in meter
%
dt = 0.01;
T = 19;
t0 = 0:dt:T;
t0_1 = t0';
n = length(t0);
m = size(t0_1);
t00 = t0;
%
fn = 0*0.05;
psi_0 = 1*0*d2r;
wn = 2*pi*fn;
wz(1) = wn;
psi(1) = psi_0 + wz(1)*t0(1);
for i = 2:n,
wz(i) = wn;
psi(i) = psi(i-1) + (wz(i) + wz(i-1))*dt/2;
end
% ====================================================
ft =1* 0.05;
wt = 2*pi*ft;
radius = (0.5*g)/(pi^2);
[x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N] = trajectory1(radius,psi_0,wt,t0);
else
    if (profile_flag ==2),
% ====================================================
% Circular motion with psi = 2*pi*fn*t
% =====================================================
% positions of two ratio sensors's positions
% ================================
xr1 = 0;% in meter
yr1 = 50;% in meter
xr2 = 0;
yr2 = -50;
h = 4;% ceiling high in meter
%
dt = 0.01;
T = 200;
t0 = 0:dt:T;
t0_1 = t0';
n = length(t0);
m = size(t0_1);
t00 = t0;
%
fn = 1*0.05;
psi_0 = 0*45*d2r;
wn = 2*pi*fn;
wz(1) = wn;
psi(1) = psi_0 + wz(1)*t0(1);
for i = 2:n,
wz(i) = wn;
psi(i) = psi(i-1) + (wz(i) + wz(i-1))*dt/2;
if (psi(i) >= 2*pi),
    psi(i) = psi(i) - 2*pi;
else
    psi(i) = psi(i);
end
end
radius = 50;
[x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N] = trajectory (radius,psi,wn);%round cycle
% ====================================================
    else
        if (profile_flag ==3),
% ===================================================
% Race track motion
% ===================================================
% positions of two ratio sensors's positions
% ================================
xr1 = 0;% in meter
yr1 = 189/2;% in meter
xr2 = 0;
yr2 = -189/2;
h = 4;% ceiling high in meter
%
v = 10;                            % ground speed in m/sec
roll_a = 1.25*d2r;                 % band-band acceleration in m/sec^2
Ts = 4 + 0.0;
T1 = 14.3 + 0.0;
[x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N,roll_angle,yaw_angle,yaw_rate,time] = trajectory2(v,roll_a,Ts,T1);
%
dt = 0.01;
T = time(length(time)-69);
t0 = 0:dt:T;
t0_1 = t0';
n = length(t0);
m = size(t0_1);
t00 = time;
psi = yaw_angle;
wz = yaw_rate;
for i = 2:n,
if (psi(i) >= 2*pi),
    psi(i) = psi(i) - 2*pi;
else
    psi(i) = psi(i);
end
end
%
        end
    end
end
% ====================================================
% Need to convert into body frame (B-frame) for accelerometer sensings and optical
% flow sensing since they are mounted on the body frame
% ====================================================
x_v_B = zeros(1,n);
y_v_B = zeros(1,n);
%
x_a_B = zeros(1,n);
y_a_B = zeros(1,n);
% ================================
for i = 1:n,
%    x_v_B(i) = [cos(psi(i)) sin(psi(i))]*[x_v_N(i);y_v_N(i)];
%    y_v_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_v_N(i);y_v_N(i)];

    x_a_B(i) = [cos(psi(i)) sin(psi(i))]*[x_a_N(i)+wz(i)*y_v_N(i);y_a_N(i)-wz(i)*x_v_N(i)];
    y_a_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_a_N(i)+wz(i)*y_v_N(i);y_a_N(i)-wz(i)*x_v_N(i)];
end
    x_v_B(1) = [cos(psi(1)) sin(psi(1))]*[x_v_N(1);y_v_N(1)];
    y_v_B(1) = [-sin(psi(1)) cos(psi(1))]*[x_v_N(1);y_v_N(1)];
for i = 2:n,
    x_v_B(i) = x_v_B(i-1) + (x_a_B(i) + x_a_B(i-1))*dt/2;
    y_v_B(i) = y_v_B(i-1) + (y_a_B(i) + y_a_B(i-1))*dt/2;
end
% ========================================================
% Define inertial sensor parameters "accelerate mpu9150"
% ========================================================
% gyroscope
if (format_1 ==1)
bz0=1*-4.37*d2r;%1*-0.89426*d2r;         % initial gyro bias in rad/sec  straight line
elseif (format_1 ==2)
bz0=1*-4.37*d2r;          % initial gyro bias in rad/sec  square
end
%=============================================================================================================
% gyroscope(bias & noise)
%=====================================
sig_arw_0 = 1*0.02;                       % gyro angle random walk = 0.02 deg/sqrt(sec)
sig_rrw_0 = 0.02/3600;                    % gyro rate random walk = 0.02 deg/3600/sec-sqrt(sec);
[bz]=Biasbg1(dt,n,m,bz0,d2r,sig_rrw_0);
[wzm]=Wzm1(wz,bz,m,n,d2r,sig_arw_0);

%=============================================================================================================
% accelerometer (biases and noises)
%=====================================
bx0=1*0.1*g;                             % initial accel bias in g in along-direction
by0=-0.1*g;                              % initial accel bias in g in perpenticular-direction
err_factor = 1.0;
%
sig_xr_0 = err_factor*0.1*g/3600;         % accel bias stability in g/sec-sqrt(sec) in along-direction
sig_yr_0 = err_factor*0.1*g/3600;         % accel bias stability in g/sec-sqrt(sec) in penpenticular-direction

% accelerate (noise)
sig_bx_0 = err_factor*0.01*g;              % accel noise in g in along-direction
sig_by_0 = err_factor*0.01*g;              % accel noise in g in penpenticular-direction

% accelerate (calculator bias)
[bx]=Biasba1(dt,n,m,bx0,sig_xr_0);
[by]=Biasbp1(dt,n,m,by0,sig_yr_0);
[axm,aym]=transform_m(x_a_B,bx,y_a_B,by,m,n,sig_bx_0,sig_by_0);

% Define ""optical flow"" parameters (noise)
opticalflow_err_factor = 1.0;
sig_x=opticalflow_err_factor*0.1;              % optical flow measurement noise in meter/sec
sig_y=opticalflow_err_factor*0.1;              % optical flow measurement noise in meter/sec in y-direction
% =========================================
% optical flow sensor model - in B-frame
% =========================================
[vxm,vym,nvx,nvy] = XYvm(x_v_B,y_v_B,m,n,sig_x,sig_y);
% calculator "Q" use bias & noise
sig_bx=sig_bx_0;%(noise)
sig_by=sig_by_0;%(noise)
sig_xr=sig_xr_0;%(bias)
sig_yr=sig_yr_0;%(bias)

% Define ""radio sensor"" parameters (noise)
radiosensor_err_factor = 1.0;
sig_x_r=radiosensor_err_factor*0.1;              % radio sensor measurement noise in meters x-direction
sig_y_r=radiosensor_err_factor*0.1;              % radio sensor measurement noise in meters y-direction
%=========================================================
[R1m,R2m,nvx_r,nvy_r] = radio_sensor_m(xr1,yr1,xr2,yr2,x_p_N,y_p_N,h,sig_x_r,sig_y_r,n,m);%4-state
%=========================================================
delta_t = dt;                                   % delta time for simulating the true dynamics = 0.01 sec
delta_s = 1*delta_t;                           % sampling at every 0.1 second for the Kalman filter
%================================================================
[sensor_step,propagation_step]=propagate_step(T,delta_t,delta_s);
%================================================================
% Define the initial conditions for the inertial sensors and the navigation
% states
%===============================================================
% Introduce initial position and velocity estimate error in N-frame
xverr = 0.1;        % in meters/sec
yverr = -0.1;       % in meters/sec
xperr = 0.5;        % in m
yperr = -0.5;       % in m
xaerr = 0;          % in m/sec^2
yaerr = 0;          % in m/sec^2
psierr = 0.5;       % in deg
[xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,by_h,bx_h,wzm_h,psi_h,bz_h]=initial_estimate_value8_radio(m,x_a_N,y_a_N,x_v_N,y_v_N,x_p_N,y_p_N,wzm,axm,aym,psi,d2r,xverr,yverr,xperr,yperr,xaerr,yaerr,psierr,y_p_N(1),x_p_N(1),bz0);
%================================================================
% Define the initial conditions for the 8-state Kalman Filter
% ==============================================================
% It is noted that we need to increase (100 times) the initial covariance matrix values
% in P00_z(8,8) to make the gyro bias error more observable
% ==============================================================
[xz_h,P00_z]=define_initial_condition_8(bx0,by0,bz0,xperr,yperr,xverr,yverr,psierr,d2r);
% ==============================================================
% for 8-state filter
F_z=zeros(8);
%F_z = [0 1 0 0 0 0 0 0
%       0 0 f23 0 f25 f26 f27 f28
%       0 0 0 0 0 0 0 0
%       0 0 0 0 1 0 0 0
%       0 f52 f53 0 0 f56 f57 f58
%       0 0 0 0 0 0 0 0
%       0 0 0 0 0 0 0 -1
%       0 0 0 0 0 0 0 0];
F_z(1,2) = 1;
F_z(4,5) = 1;
F_z(7,8) = -1;
Q_z=zeros(8);

% ============================================================
% Start the simulation run
% ============================================================
k=2;
for i=1:sensor_step
	for j=1:propagation_step
        k=1+j+(i-1)*propagation_step;
        
        bx_h(k)=bx_h(k-1);
        by_h(k)=by_h(k-1);
        bz_h(k)=bz_h(k-1);
% =========================================        
% Perform inertial navigation computations
% =========================================
        [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,wzm_h,psi_h]=inertial_navigation_computation8_radio(xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,xam_Nh,yam_Nh,axm,aym,bx_h,by_h,psi_h,wzm,bz_h,k,dt,0);
    
% ===========================================================
% Perform Kalman filter propagation for 8-state Kalman filter
% ===========================================================
        [phi_z,Q_z,F_z]=define_Dymamic_equation8_radio(F_z,Q_z,axm_h,aym_h,xvm_Nh,yvm_Nh,wzm_h,psi_h,sig_bx,sig_by,sig_xr,sig_yr,sig_arw_0,sig_rrw_0,dt,k,1);
        [xz_h,P00_z]=Kalman_Filter_estimate1_radio(xz_h,phi_z,P00_z,Q_z,dt);
    end                 % end of filter propagation step
% =======================================================
% Perform Kalman filter updates for 8-state filter
% =======================================================
    [H,R,R1m_h,R2m_h]=radio_discrete_8_2_EKF(xr1,yr1,xr2,yr2,xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,h,k);
%   
    [P00_z,K_z,z_update]=Kalman_Filter_update_8_2_radio(P00_z,R1m,R2m,H,R,R1m_h,R2m_h,k);
    
    [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,psi_h,bz_h,bx_h,by_h]=upon_radiosensor_measurement_8_2(xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,k,psi_h,bz_h,bx_h,by_h,z_update);
    
% reset the errors after the filter updates
    xz_h=zeros(8,1);
end                     % end of one Monte Caro run
% ==================================================================================
% Store the filter performance data: accelerometers' bias estimate errors and gyro
% bias estimate error
% ==================================================================================
m1 = 12500;
m2 = k-1;
bxerrave = mean((bx(m1:m2)-bx_h(m1:m2))/g);
byerrave = mean((by(m1:m2)-by_h(m1:m2))/g);
bzerrave = mean((bz(m1:m2)-bz_h(m1:m2))*r2d);
bxerrstd = std((bx(m1:m2)-bx_h(m1:m2))/g);
byerrstd = std((by(m1:m2)-by_h(m1:m2))/g);
bzerrstd = std((bz(m1:m2)-bz_h(m1:m2))*r2d);
bxerr(ii) = abs(bxerrave) + 3*bxerrstd;
byerr(ii) = abs(byerrave) + 3*byerrstd;
bzerr(ii) = abs(bzerrave) + 3*bzerrstd;
%
xperrave(ii) = abs(mean(x_p_N(m1:m2)-xpm_Nh(m1:m2)'));
yperrave(ii) = abs(mean(y_p_N(m1:m2)-ypm_Nh(m1:m2)'));
xperrstd(ii) = std(x_p_N(m1:m2)-xpm_Nh(m1:m2)');
yperrstd(ii) = std(y_p_N(m1:m2)-ypm_Nh(m1:m2)');

%xperr(ii) = xperrave(ii) + 3*xperrstd(ii),
%yperr(ii) = yperrave(ii) + 3*yperrstd(ii),
end                     % end of N Monte Caro runs
%
for k1 = 1:ii,
    xperr(k1) = xperrave(k1) + 3*xperrstd(k1);
    yperr(k1) = yperrave(k1) + 3*yperrstd(k1);
end
%
bxerrp = mean(bxerr);
byerrp = mean(byerr);
bzerrp = mean(bzerr);
xperrp = mean(xperr);
yperrp = mean(yperr);
[bxerrp byerrp bzerrp],
[xperrp yperrp],

% ============================================================
% Plot the performance
% ============================================================
% plot the filter performance related to the 2-state Kalman filters
%t00 = time;
plot18;
% ============================================================
% plot the filter performance related to the 4-state Kalman filter
plot28;
%
