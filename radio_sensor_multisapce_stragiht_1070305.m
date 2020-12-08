% Program name: radio_sensor_8_4_v1.m
% A new program that deal with non-constant psi angle where
% a 8-state Kalman filter with 4 radio sensors without optical flow
% sensors
close all;
clear all;
clc;

dwm_trans_10anchors
% dwm_trans_3

format short e       
%format short e
r2d = (180/pi);
d2r = (pi/180);
g =9.8;
scale_factor_err = -0.0512;%; °f -0.0512 ¶¶ -0.2577-0.06112
format_1 = 1; % stiaight line = 1 (-0.89426), square = 2 (-4.37)  
i_1=1;i_2=1;i_3=1;i_4=1;i_5=1;
i_21=1;i_22=1;i_23=1;i_24=1;i_25=1;i_26=1;
%

for ii = 1:100,
%======================================================
% Design parameters to be tuned
%=====================================================
% trajectory is generated in navigation frame (N-frame)
%================================
% Set motion profile flags

profile_flag = 1;
%
if (profile_flag ==1),
% =====================================================
%(y_anchor_2 = 9.37; x_anchor_2 = 0.55;z_anchor_2 = 1.83;)
%unit in meter
%             idx x y z 
coordinate = [1,3.65,-0.705,2.515;
              2,3.672,-6.565,2.46;
              3,11.144,-0.705,2.44;
              4,11.147,-6.565,2.86;
              5,22.428,-0.705,2.55;
              6,22.422,-6.565,2.47;
              7,0,0.9,2.75;
              8,0.17,11.8,2.72;
              9,7.78,11.7,2.5;
              10,7.65,0.13,2.7];
%h_i is ceiling high 
xr1 = coordinate(7,2);
yr1 = coordinate(7,3);
h_1 = coordinate(7,4); 

xr2 = coordinate(8,2);   
yr2 = coordinate(8,3); 
h_2 = coordinate(8,4); 

xr3 = coordinate(9,2);% in meter7.85¡@¡@7.78 7.85
yr3 = coordinate(9,3);
h_3 = coordinate(9,4);

xr4 = coordinate(10,2);
yr4 = coordinate(10,3);
h_4 = coordinate(10,4);

dt = 0.01;
T = 39.9;
t0 = 0:dt:T;
t0_1 = t0';
n = length(t0);
m = size(t0_1);
t00 = t0;
%
fn = 0*0.05;
psi_0 = 1*0*d2r;
wn = 0;
wz(1) = wn;
psi(1) = psi_0 + wz(1)*t0(1);
for i = 2:n,
wz(i) = wn;
psi(i) = psi(i-1) + (wz(i) + wz(i-1))*dt/2;
end
% ====================================================
ft =1* 0.5;
wt = 2*pi*ft;
% radius = 25;
radius = (0.5*g)/(pi^2);
[x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N] = trajectory1(radius,psi_0,wt,t0);

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
% [wzm]=Wzm1(wz,bz,m,n,d2r,sig_arw_0);
% k=2
% for j=1:propagation_step
% psi_(k) = psi_(k-1) + wzm(k)*dt;
% if (psi_(k)>= 2*pi),
%     psi_(k) = psi_(k) - 2*pi;
% else
%     psi_(k) = psi_(k);
% end
% xam_(k) = cos(psi(k-1))*axm_h(k-1)-sin(psi(k-1))*aym(k-1)-wzm(k-1)*yvm_Nh(k-1);
% yam_(k) = sin(psi(k-1))*axm_h(k-1)+cos(psi(k-1))*aym(k-1)+wzm(k-1)*xvm_Nh(k-1);
% xvm_(k) = xvm(k-1)+(xam_(k)+xam_(k-1))*dt/2;
% yvm_(k) = yvm(k-1)+(yam_(k)+yam_(k-1))*dt/2;
% end
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
sig_bx_0 = err_factor*0.1*g;              % accel noise in g in along-direction
sig_by_0 = err_factor*0.1*g;              % accel noise in g in penpenticular-direction

% accelerate (calculator bias)
[bx]=Biasba1(dt,n,m,bx0,sig_xr_0);
[by]=Biasbp1(dt,n,m,by0,sig_yr_0);
[axm,aym]=transform_m(x_a_B,bx,y_a_B,by,m,n,sig_bx_0,sig_by_0);

% calculator "Q" use bias & noise
sig_bx=sig_bx_0;%(noise)
sig_by=sig_by_0;%(noise)
sig_xr=sig_xr_0;%(bias)
sig_yr=sig_yr_0;%(bias)

% Define ""4 radio sensor"" parameters (noise)
radiosensor_err_factor = 1.0;
sig_x_r=radiosensor_err_factor*0.1;              % radio sensor measurement noise in meters x-direction
sig_y_r=radiosensor_err_factor*0.1;              % radio sensor measurement noise in meters y-direction
%=========================================================
% [R1m,R2m,R3m,R4m,nvx_r,nvy_r] = radio_sensor_m_4(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,x_p_N,y_p_N,h_1,h_2,h_3,h_4,sig_x_r,sig_y_r,n,m);%4 four radio sensors
% [Rm_data,Rm_Index,nvx_r,nvy_r] = radio_sensor_m_mult(coordinate,x_p_N,y_p_N,sig_x_r,sig_y_r,n,m);%8-state
%=========================================================
delta_t = dt;                                   % delta time for simulating the true dynamics = 0.01 sec
delta_s = 2*delta_t;                           % sampling at every 0.1 second for the Kalman filter
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
psierr = 1;       % in deg
[xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,by_h,bx_h,wzm_h,psi_h,bz_h]=initial_estimate_value8_radio(m,x_a_N,y_a_N,x_v_N,y_v_N,x_p_N,y_p_N,gro(3,:)',acc(1,:)',acc(2,:)',psi,d2r,xverr,yverr,xperr,yperr,xaerr,yaerr,psierr,position(1,1),position(2,1),bz0);%gro ,acc(1,:)',acc(2,:)'
% [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,by_h,bx_h,wzm_h,psi_h,bz_h]=initial_estimate_value8_radio(m,x_a_N,y_a_N,x_v_N,y_v_N,x_p_N,y_p_N,wzm,low_gro(3,:)',low_acc(1,:)',low_acc(2,:)',d2r,xverr,yverr,xperr,yperr,xaerr,yaerr,psierr,position(1,1),position(2,1));%gro ,acc(1,:)',acc(2,:)'

%================================================================
% Define the initial conditions for the 8-state Kalman Filter
% ==============================================================
% It is noted that we increase (100 times) the initial covariance matrix values
% in P00_z(8,8) to make the gyro bias error more observable, and in
% P00_z(3,3), P00_z(6,6) to make accel bias errors converge faster
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
l=1;
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
%         [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,wzm_h,psi_h]=inertial_navigation_computation8_radio(xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,xam_Nh,yam_Nh,low_acc(1,:)',low_acc(2,:)',bx_h,by_h,psi_h,low_gro(3,:)',bz_h,k,dt,scale_factor_err);
        [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,xam_Nh,yam_Nh,axm_h,aym_h,wzm_h,psi_h]=inertial_navigation_computation8_radio(xvm_Nh,yvm_Nh,xpm_Nh,ypm_Nh,xam_Nh,yam_Nh,acc(1,:)',acc(2,:)',bx_h,by_h,psi_h,gro(3,:)',bz_h,k,dt,scale_factor_err);%gro(3,:)' acc(1,:)',acc(2,:)'
    
% ===========================================================
% Perform Kalman filter propagation for 8-state Kalman filter
% ===========================================================
    
        [phi_z,Q_z,F_z]=define_Dymamic_equation8_radio(F_z,Q_z,axm_h,aym_h,xvm_Nh,yvm_Nh,wzm_h,psi_h,sig_bx,sig_by,sig_xr,sig_yr,sig_arw_0,sig_rrw_0,dt,k,format_1);
   
        [xz_h,P00_z]=Kalman_Filter_estimate1_radio(xz_h,phi_z,P00_z,Q_z,dt);
    end                 % end of filter propagation step
% =======================================================
% Perform Kalman filter updates for 8-state filter
% =======================================================
   
    [H,R,R1m_h,R2m_h,R3m_h,R4m_h]=radio_discrete_8_4_EKF(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,h_1,h_2,h_3,h_4,k,format_1);
   
%     if (k<2001)
%     c=2;%,Mu_count,i_1,i_2,i_3,i_4,i_5,i_1,i_2,i_3,i_4,i_5
%     elseif((2000< k) && (k <2500))
%     c=0.01;%,Mu_count,i_1,i_2,i_3,i_4,i_5,i_1,i_2,i_3,i_4,i_5
%     else
%     c=0.01;
%     end
    
    
    c=0.5;%,Mu_count,i_1,i_2,i_3,i_4,i_5,i_1,i_2,i_3,i_4,i_5
%     [P00_z,K_z,z_update,zxm_z1,zxm_z2,zxm_z3,zxm_z4,Mu_count_5,Mu_count_4,Mu_count_3,Mu_count_2,Mu_count_1,Mu_count_21,Mu_count_22,Mu_count_23,Mu_count_24,Mu_count_25,Mu_count_26,Mu_count_31,Mu_count_32,Mu_count_33,Mu_count_34,i_1,i_2,i_3,i_4,i_5,i_21,i_22,i_23,i_24,i_25,i_26]=Kalman_Filter_update_8_4_radio_v1(P00_z,pos(1,:)',pos(2,:)',pos(3,:)',pos(4,:)',H,R,R1m_h,R2m_h,R3m_h,R4m_h,k,c,i_1,i_2,i_3,i_4,i_5,i_21,i_22,i_23,i_24,i_25,i_26);

    [P00_z,K_z,z_update,zxm_z1,zxm_z2,zxm_z3,zxm_z4]=Kalman_Filter_update_8_4_radio_v1(P00_z,pos(1,:)',pos(2,:)',pos(3,:)',pos(4,:)',H,R,R1m_h,R2m_h,R3m_h,R4m_h,k,c);
%     [P00_z,K_z,z_update]=Kalman_Filter_update_8_4_radio(P00_z,pos(1,:)',pos(2,:)',pos(3,:)',pos(4,:)',H,R,R1m_h,R2m_h,R3m_h,R4m_h,k,c);
    [xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,psi_h,bz_h,bx_h,by_h]=upon_radiosensor_measurement_8_2(xpm_Nh,ypm_Nh,xvm_Nh,yvm_Nh,k,psi_h,bz_h,bx_h,by_h,z_update);
    
% reset the errors after the filter updates
    xz_h=zeros(8,1);
    Mu_z_1(l,:) = zxm_z1(1,k);
    Mu_z_2(l,:) = zxm_z2(1,k);
    Mu_z_3(l,:) = zxm_z3(1,k);
    Mu_z_4(l,:) = zxm_z4(1,k);
    
%     Mu_count__1(1,l) = Mu_count_1;
%     Mu_count__1(2,l) = Mu_count_2;
%     Mu_count__1(3,l) = Mu_count_3;
%     Mu_count__1(4,l) = Mu_count_4;
%     Mu_count__1(5,l) = Mu_count_5;
%     
%     Mu_count__2(1,l) = Mu_count_21; % not anchor1 anchor2
%     Mu_count__2(2,l) = Mu_count_22; % not anchor1 anchor3
%     Mu_count__2(3,l) = Mu_count_23; % not anchor1 anchor4
%     Mu_count__2(4,l) = Mu_count_24; % not anchor2 anchor3
%     Mu_count__2(5,l) = Mu_count_25; % not anchor2 anchor4
%     Mu_count__2(6,l) = Mu_count_26; % not anchor3 anchor4
%     
%     Mu_count__3(1,l) = Mu_count_31; % not % Not Anchor1 Anchor2 Anchor3
%     Mu_count__3(2,l) = Mu_count_32; % not % Not Anchor1 Anchor2 Anchor4
%     Mu_count__3(3,l) = Mu_count_33; % not % Not Anchor1 Anchor3 Anchor4
%     Mu_count__3(4,l) = Mu_count_34; % not % Not Anchor2 Anchor3 Anchor4
   
    l = l+1;
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
plot48;
% ============================================================
% plot the filter performance related to the 4-state Kalman filter
plot28;
%
