% Program name: radio_sensor_2.m
% A new program that deal with non-constant psi angle where
% a 6-state Kalman filter is needed
close all;
clear all;
clc;

r2d = (180/pi);
d2r = (pi/180);
g =9.8;
%
for ii = 1:100,
dt = 0.01;% delta time for simulating the true dynamics
T=200;
t0 = 0:dt:T;
t0_1 = t0';
n = length(t0);
m = size(t0_1);
%======================================================
% Design parameters to be tuned
%=====================================================
% (1) positions of two ratio sensors's positions
%================================
xr1 = 0;% in meter
yr1 = 50;% in meter
xr2 = 0;
yr2 = -50;
h = 4;% ceiling high in meter
%================================
% (2) motion profiles
% trajectory is generated in navigation frame (N-frame)
%================================ 
fn = 0*0.05;
psi_0 = 1*45*d2r;
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
%
% ====================================================
% (a) Circular motion with psi = 2*pi*fn*t
%=====================================================
%radius = 50;
%[x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N] = trajectory (radius,psi,wn);%round cycle
% ====================================================
% (b) Straight motion with psi = constant angle, or wn = 0
% ====================================================
ft =1* 0.05;
wt = 2*pi*ft;
radius = 50;
[x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N] = trajectory1(radius,psi_0,wt,t0);
% ====================================================
% Need to convert into body frame (B-frame) for accelerometer sensings and optical
% flow sensing since they are mounted on the body frame
% From N-frame to B-frame
%[x_p_B,x_v_B,x_a_B,y_p_B,y_v_B,y_a_B] = NtoB_transfer(x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N,psi,n);
% ================================================================================================

x_v_B = zeros(1,n);
y_v_B = zeros(1,n);

x_a_B = zeros(1,n);
y_a_B = zeros(1,n);
% ================================
% it in noted that the following velocity transformation is true only if
% psi angle is constant, or wn is equal to zero
% It is also noted that for the current N-frame trajectory motion: circular
% motion, the following velocity transformation is correct - this is a
% special case
% The body accelerations include the Coriolis effect force: w x V_N
% ================================
for i = 1:n
    x_v_B(i) = [cos(psi(i)) sin(psi(i))]*[x_v_N(i);y_v_N(i)];
    y_v_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_v_N(i);y_v_N(i)];

    x_a_B(i) = [cos(psi(i)) sin(psi(i))]*[x_a_N(i)+wz(i)*y_v_N(i);y_a_N(i)-wz(i)*x_v_N(i)];
    y_a_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_a_N(i)+wz(i)*y_v_N(i);y_a_N(i)-wz(i)*x_v_N(i)];
end
% ========================================================
% Define inertial sensor parameters "accelerate mpu9150"
% ========================================================
% gyroscope
bz0=1*0.5*d2r;                            % initial gyro bias in rad/sec
%=============================================================================================================
% gyroscope(bias & noise) %4-state
%=====================================
sig_arw_0 = 1*0.02;                       % gyro angle random walk = 0.02 deg/sqrt(sec)
sig_rrw_0 = 0.02/3600;                    % gyro rate random walk = 0.02 deg/3600/sec-sqrt(sec);
[bz]=Biasbg1(dt,n,m,bz0,d2r,sig_rrw_0);
[wzm]=Wzm1(wz,bz,m,n,d2r,sig_arw_0); %%%%%%%%%%%%%%%%% not value

%=============================================================================================================
% accelerometer (biases and noises)
%=====================================
ba0=1*0.1*g;                             % initial accel bias in g in along-direction
bp0=-0.1*g;                              % initial accel bias in g in perpenticular-direction
err_factor = 1.0;
%
sig_xr_0 = err_factor*0.1*g/3600;         % accel bias stability in g/sec-sqrt(sec) in along-direction
sig_yr_0 = err_factor*0.2*g/3600;         % accel bias stability in g/sec-sqrt(sec) in penpenticular-direction

% accelerate (noise)
sig_bx_0 = err_factor*0.01*g;              % accel noise in g in along-direction
sig_by_0 = err_factor*0.02*g;              % accel noise in g in penpenticular-direction

% accelerate (calculator bias)
[bx]=Biasba1(dt,n,m,ba0,sig_xr_0);
[by]=Biasbp1(dt,n,m,bp0,sig_yr_0);% 到t秒時的bias
[axm,aym]=transform_m(x_a_B,bx,y_a_B,by,m,n,sig_bx_0,sig_by_0);% 利用 bias & noise 得到量測的 aam & apm

% Define ""optical flow"" parameters (noise)
opticalflow_err_factor = 1.0;
% % % % % % % % % % sig_x=opticalflow_err_factor*0.01;              % optical flow measurement noise in meter/sec
% % % % % % % % % % sig_y=opticalflow_err_factor*0.02;              % optical flow measurement noise in meter/sec in y-direction
%
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

% Define ""radio sensor"" parameters (noise) %4-state
radiosensor_err_factor = 1.0;
sig_x_r=radiosensor_err_factor*0.1;              % radio sensor measurement noise in meters x-direction
sig_y_r=radiosensor_err_factor*0.1;              % radio sensor measurement noise in meters y-direction
%=========================================================
% [Px_N,Py_N] = BtoN_transfer(x_p,y_p,x_v,y_v);%transfer
%[x_v_N,y_v_N,x_p_N,y_p_N] = BtoN_transfer(x_v,y_v,x_p,y_p,psi,n);%4-state
% It is noted that we don't need the above BtoN_transfer function since the
% x_p_N, y_p_N have been generated from the function: trajector.m. In
% addition, the above transfermation in general is not correct unless the
% psi angle is constant
%=========================================================
[R1m,R2m,nvx_r,nvy_r] = radio_sensor_m(xr1,yr1,xr2,yr2,x_p_N,y_p_N,h,sig_x_r,sig_y_r,n,m);%4-state
%=========================================================
delta_t = dt;                           % delta time for simulating the true dynamics = 0.01 sec
%delta_s = delta_t;                      % sampling at every 0.01 second for the Kalman filter
delta_s = 10*delta_t;                      % sampling at every 0.01 second for the Kalman filter
%================================================================
[sensor_step,propagation_step]=propagate_step(T,delta_t,delta_s);
%================================================================
% Define the initial conditions for the inertial sensors
% one for the two 2-state kalman filters, and the other for the 4-state kalman
% filter
%================================================================
% Introduce +- 0.1m/sec initial errors in Vx_h and Vy_h; 0.1 deg error in
% psi_h
Vxerr = 0.1;
Vyerr = -0.1;
%===============================================================
% 2-state kalman filters
[axm_h,aym_h,Vx_h,Vy_h,vxm_h,vym_h,bx_h,by_h]=initial_estimate_value1(m,vxm,vym,axm,aym,x_v_B,y_v_B,Vxerr,Vyerr);
%[axm_h,aym_h,Vx_h,Vy_h,vxm_h,vym_h,bx_h,by_h]=initial_estimate_value1(m,vxm,vym,axm,aym,x_v_B,y_v_B,Vxerr,Vyerr);
% 4-state kalman filter
% Introduce initial position and velocity estimate error in N-frame
xverr = 0.1;        % in meters
yverr = -0.1;       % in meters
xperr = 0.5;        % in m/sec
yperr = -0.5;       % in m/sec
xaerr = 0;          % in m/sec^2
yaerr = 0;          % in m/sec^2
psierr = 0.5;       % in deg
[xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,xam_Nh,yam_Nh,wzm_h,psi_h,bz_h]=initial_estimate_value2_radio(m,x_a_N,y_a_N,x_v_N,y_v_N,x_p_N,y_p_N,wzm,psi,d2r,xverr,yverr,xperr,yperr,xaerr,yaerr,psierr);
%================================================================
% Define the initial conditions for the Kalman Filter
% one for the two 2-state kalman filters, and the other for the 6-state kalman
% filter
%===============================================================
[x_h,y_h,P00_x,P00_y]=define_initial_condition1(Vxerr,Vyerr,ba0,bp0);
%[x_h,y_h,P00_x,P00_y]=define_initial_condition1(Vxerr,Vyerr,ba0,bp0);
% ==============================================================
% It is noted that we need to increase (100 times) the initial covariance matrix values
% in P00_z(6,6) to make the gyro bias error more observable
% ==============================================================
[xz_h,P00_z]=define_initial_condition_2(bz0,xperr,yperr,xverr,yverr,psierr,d2r);
% ==============================================================
% for 2-state filter
H_x = zeros(1,2);H_x(1,1) = 1;
H_y = zeros(1,2);H_y(1,1) = 1;
R_x = sig_x^2;
R_y = sig_y^2;
F = [0 -1
     0  0];
% for 4-state filter
F_z=zeros(6);
%F_z = [0 1 0 0 0 0
%       0 0 0 f24 f25 f26
%       0 0 0 1 0 0
%       0 f42 0 0 f45 f46
%       0 0 0 0 0 -1
%       0 0 0 0 0 0];
F_z(1,2) = 1;
F_z(3,4) = 1;
F_z(5,6) = -1;
Q_z=zeros(6);

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
        %
% for two 2-state Kalman filters
        [axm_h,aym_h,Vx_h,Vy_h]=inertial_navigation_computation1(axm_h,aym_h,Vx_h,Vy_h,axm,aym,bx_h,by_h,k,dt);
%        [axm_h,aym_h,Vx_h,Vy_h]=inertial_navigation_computation1(axm_h,aym_h,Vx_h,Vy_h,axm,aym,bx_h,by_h,k,dt);
% for 4-state Kalman filter where Vx_h and Vy_h produced from
% 2-state Kalman filters are used in the 4-state Kalman filter
% ==============================================================
% for debuging purpose we subtitute the true body velocities for
% the estimated velocities from the two 2-state Kalman filter
% xvm_Nh(k) = cos(psi_h(k))*Vx_h(k)-sin(psi_h(k))*Vy_h(k);
% yvm_Nh(k) = sin(psi_h(k))*Vx_h(k)+cos(psi_h(k))*Vy_h(k);
% xvm_Nh(k) = cos(psi_h(k))*x_v_B(k)-sin(psi_h(k))*y_v_B(k);
% yvm_Nh(k) = sin(psi_h(k))*x_v_B(k)+cos(psi_h(k))*y_v_B(k);
        %
        [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,wzm_h,psi_h]=inertial_navigation_computation2_radio(xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,wzm_h,psi_h,wzm,bz_h,k,dt);
    
% for 2-state Kalman filter
        [phi,Q_x,Q_y]=define_Dymamic_equation1(F,sig_bx,sig_by,sig_xr,sig_yr,dt);
% for 6-state Kalman filter
        [phi_z,Q_z,F_z]=define_Dymamic_equation2_radio(F_z,Q_z,axm_h,aym_h,xvm_Nh,yvm_Nh,wzm_h,psi_h,P00_x,P00_y,sig_arw_0,sig_rrw_0,dt,k);
      
        [x_h,y_h,P00_x,P00_y]=Kalman_Filter_estimate1(x_h,y_h,phi,P00_x,P00_y,Q_x,Q_y,dt);%update
        %[x_h,y_h,P00_x,P00_y]=Kalman_Filter_estimate1(x_h,y_h,phi,P00_x,P00_y,Q_x,Q_y,dt);
        %
        [xz_h,P00_z]=Kalman_Filter_estimate1_radio(xz_h,phi_z,P00_z,Q_z,dt);%6-state
    end                 % end of filter propagation step
    

    % for 2-state Kalman filter update
    [x_h,y_h,P00_x,P00_y,K_x,K_y,x_update,y_update]=Kalman_Filter_update1(x_h,y_h,P00_x,P00_y,Vx_h,Vy_h,vxm,vym,H_x,H_y,R_x,R_y,k);%update
    %[x_h,y_h,P00_x,P00_y,K_x,K_y,x_update,y_update]=Kalman_Filter_update1(x_h,y_h,P00_x,P00_y,Vx_h,Vy_h,vxm,vym,H_x,H_y,R_x,R_y,k);
    % for 4-state Kalman filter update
    %H,F_z,
    [H,R,R1m_h,R2m_h]=radio_discrete_2_EKF(xr1,yr1,xr2,yr2,xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,h,k);%6-state
%   
    [P00_z,K_z,z_update]=Kalman_Filter_update2_radio(xz_h,P00_z,R1m,R2m,H,R,R1m_h,R2m_h,k);%6-state

    [Vx_h,Vy_h,bx_h,by_h]=upon_optimalflow_measurement(x_update,y_update,k,Vx_h,Vy_h,bx_h,by_h);
    %[Vx_h,Vy_h,bx_h,by_h]=upon_optimalflow_measurement(x_update,y_update,k,Vx_h,Vy_h,bx_h,by_h);
    
    [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,psi_h,bz_h]=upon_radiosensor_measurement_2(xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,k,psi_h,bz_h,z_update);%6-state
    
% reset the errors after the filter updates
    x_h=zeros(2,1);
    y_h=zeros(2,1);
    xz_h=zeros(6,1);%4-state
end                     % end of one Monte Caro run
% Store the filter performance data: accelerometers' bias estimate errors and gyro
% bias estimate error
m1 = 10000;
m2 = k-1;
bxerrave = mean((bx(m1:m2)-bx_h(m1:m2))/g);
byerrave = mean((by(m1:m2)-by_h(m1:m2))/g);
bzerrave = mean((bz(m1:m2)-bz_h(m1:m2))*r2d);
bxerrstd = std((bx(m1:m2)-bx_h(m1:m2))/g);
byerrstd = std((by(m1:m2)-by_h(m1:m2))/g);
bzerrstd = std((bz(m1:m2)-bz_h(m1:m2))*r2d);
bxerr(ii) = bxerrave + 3*bxerrstd;
byerr(ii) = byerrave + 3*byerrstd;
bzerr(ii) = bzerrave + 3*bzerrstd;
end                     % end of N Monte Caro runs
bxerrp = mean(bxerr),
byerrp = mean(byerr),
bzerrp = mean(bzerr),

% ============================================================
% Plot the performance
% ============================================================
% plot the filter performance related to the 2-state Kalman filters
plot1;
% ============================================================
% plot the filter performance related to the 4-state Kalman filter
plot2;
%
