% Program name:RTKandIMU_positioning.m
% A new program that deal with non-constant psi angle where
% a 8-state Kalman filter with GPS sensor

close all;
clear all;
clc;

% import real data
exp201910161116;

format short e
r2d = (180/pi);
d2r = (pi/180);
g =9.8;
% Set motion profile flags
profile_flag = 1;


%% Straight motion
%
if (profile_flag ==1)
    dt = 0.1;
    T = 700;
    t0 = 0:dt:T;
    t0_1 = t0';
    n = length(t0);
    m = size(t0_1);
    t00 = t0;
    %
    fn = 0*0.05;
    psi_0 = 1*45*d2r; % Slope of trajectory
    wn = 2*pi*fn;
    wz(1) = wn;
    psi(1) = psi_0 + wz(1)*t0(1);
    for i = 2:n,
    wz(i) = wn;
    psi(i) = psi(i-1)+ (wz(i) + wz(i-1))*dt/2;
    end
    % ====================================================
    ft =1* 0.05;
    wt = 2*pi*ft;
    radius_x = 50;
    radius_y = 25;
    [x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N] = trajectory_mult(radius_x,radius_y,psi_0,wt,t0);
else if (profile_flag ==2)
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
    end
end
% ====================================================
% Need to convert into body frame (B-frame) for accelerometer sensings
% since they are mounted on the body frame
% ====================================================
x_v_B = zeros(1,n);
y_v_B = zeros(1,n);
%
x_a_B = zeros(1,n);
y_a_B = zeros(1,n);
% ================================
for i = 1:n,
%     x_a_B(i) = [cos(psi(i)) sin(psi(i))]*[x_a_N(i)+wz(i)*y_v_N(i);y_a_N(i)-wz(i)*x_v_N(i)];
%     y_a_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_a_N(i)+wz(i)*y_v_N(i);y_a_N(i)-wz(i)*x_v_N(i)];
    x_a_B(i) = [cos(psi(i)) sin(psi(i))]*[x_a_N(i);y_a_N(i)];
    y_a_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_a_N(i);y_a_N(i)];
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
bz0=1*0.5*d2r;                            % initial gyro bias in rad/sec
%=============================================================================================================
% gyroscope(bias & noise)
%=====================================
sig_arw_0 = 1*0.02;                       % gyro angle random walk = 0.02 deg/sqrt(sec)
sig_rrw_0 = 1*0.02/3600;                    % gyro rate random walk = 0.02 deg/3600/sec-sqrt(sec);
[bz]=Biasbg1(dt,n,m,bz0,d2r,sig_rrw_0);
[wzm]=Wzm1(wz,bz,m,n,d2r,sig_arw_0);

%=============================================================================================================
% accelerometer (biases and noises)
%=====================================
bx0=1*0.1*g;                             % initial accel bias in g in along-direction
by0=-1*0.1*g;                              % initial accel bias in g in perpenticular-direction
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

% calculator "Q" use bias & noise
sig_bx=sig_bx_0;%(noise)
sig_by=sig_by_0;%(noise)
sig_xr=sig_xr_0;%(bias)
sig_yr=sig_yr_0;%(bias)

% Define ""radio sensor"" parameters (noise)
RTK_err_factor = 1.0;
sig_x_r=RTK_err_factor*0.05;              % RTK measurement noise in meters x-direction
sig_y_r=RTK_err_factor*0.05;              % RTK measurement noise in meters y-direction

nvx_r=normrnd(0,sig_x_r,1,n);
nvy_r=normrnd(0,sig_y_r,1,n);
Pm_data = [p_N(:,1) p_N(:,2)];
Pm_IMU_h = zeros(row,2);
% RTK sensor to back wheel center in B-frame
delta_L_B = [0 0.2];
Vehicle_B_N = zeros(2);
%=========================================================
delta_t = dt;                                   % delta time for simulating the true dynamics = 0.01 sec
delta_s = 1*delta_t;                            % sampling at every 0.1 second for the Kalman filter
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
psierr = 0.5;       % in deg                                                                                                                                                                                           %    x         y                        
[xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,by_h,bx_h,wzm_h,psi_h,bz_h]=initial_estimate_value8_radio(m,x_a_N,y_a_N,x_v_N,y_v_N,x_p_N,y_p_N,groz,accx,accy,psi,d2r,xverr,yverr,xperr,yperr,xaerr,yaerr,psierr,p_N(1,1),p_N(1,2),bz0);
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
        [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,wzm_h,psi_h]=inertial_navigation_computation8_radio(xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,xam_Nh,yam_Nh,accx,accy,bx_h,by_h,psi_h,groz,bz_h,k,dt,0);
    
% ===========================================================
% Perform Kalman filter propagation for 8-state Kalman filter
% ===========================================================
        [phi_z,Q_z,F_z]=define_Dymamic_equation8_radio(F_z,Q_z,axm_h,aym_h,xvm_Nh,yvm_Nh,wzm_h,psi_h,sig_bx,sig_by,sig_xr,sig_yr,sig_arw_0,sig_rrw_0,dt,k,1);
        [xz_h,P00_z]=Kalman_Filter_estimate1_radio(xz_h,phi_z,P00_z,Q_z,dt);
    end                 % end of filter propagation step
% =======================================================
% Perform Kalman filter updates for 8-state filter
% =======================================================

    [H,R,Pm_h]=RTK_discrete_EKF(xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,k,1);
    
    if 1
        % RTK center to back wheel center   
%         Vehicle_B_N(1,1) = cos(psi_h(k-1));
%         Vehicle_B_N(1,2) = -sin(psi_h(k-1));
%         Vehicle_B_N(2,1) = sin(psi_h(k-1));
%         Vehicle_B_N(2,2) = cos(psi_h(k-1));
%         if k<=5
            ang = attz(k);
%         else
%             ang = mean(attz(k-5:k+5));
%         end
        Vehicle_B_N(1,1) = cos(ang);
        Vehicle_B_N(1,2) = -sin(ang);
        Vehicle_B_N(2,1) = sin(ang);
        Vehicle_B_N(2,2) = cos(ang);        
        Pm_IMU_h(k,:) = Pm_data(k,:) - (Vehicle_B_N * delta_L_B')';
    else 
        Pm_IMU_h(k,:) = Pm_data(k,:);
    end
%   
    [P00_z,K_z,z_update]=Kalman_Filter_update_RTK(P00_z,Pm_IMU_h,H,R,Pm_h,k);
    
    [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,psi_h,bz_h,bx_h,by_h]=upon_radiosensor_measurement_8_2(xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,k,psi_h,bz_h,bx_h,by_h,z_update);
    
% reset the errors after the filter updates
    xz_h=zeros(8,1);
end                     % end of one Monte Caro run

% position error
delta_P = [Pm_IMU_h(1:k,1) - xpm_Nh(1:k) Pm_IMU_h(1:k,2) - ypm_Nh(1:k)];
delta_euler = psi_h(100:k) - attz(100:k); 
delta_P_mean = mean(delta_P);
delta_P_var = var(delta_P);
delta_euler_mean = mean(delta_euler)*r2d;
delta_euler_var = var(delta_euler);

% plot
% plotRTKsimul;
plotRTKexp;

