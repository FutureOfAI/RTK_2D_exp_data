close all;
clear all;
clc;

r2d = (180/pi);
d2r = (pi/180);

for ii = 1:1,

dt = 0.01;% delta time for simulating the true dynamics
T=80;
t0 = 0:dt:T;
t0_1 = t0';
% initial trajectory value
radius = 50; 
fn = 0.05;
g =9.8;
wn = 2*pi*fn;
n = length(t0);
m = size(t0_1);
for i = 1:n,
wz(i) = wn;
psi(i) = wz(i)*t0(i);
if (psi(i) >= 2*pi),
    psi(i) = psi(i) - 2*pi;
else
    psi(i) = psi(i);
end
end
%
% ====================================================
% trajectory is generated in navigation frame (N-frame)
% ====================================================
[x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N] = trajectory (radius,psi,wn);%round cycle
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
% The body accelerations include the Coriolis effect force: w x V_B
% ================================
for i = 1:n
    x_v_B(i) = [cos(psi(i)) sin(psi(i))]*[x_v_N(i);y_v_N(i)];
    y_v_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_v_N(i);y_v_N(i)];

    x_a_B(i) = [cos(psi(i)) sin(psi(i))]*[x_a_N(i);y_a_N(i)] + radius*wn^2;
    y_a_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_a_N(i);y_a_N(i)] + 0;
end
% ======================================================================
%[x_p,x_v,x_a,y_p,y_v,y_a] = trajectory (radius,wn,t0);%round cycle 
%figure (9)
%subplot(211)
%plot(t0,x_a_B/g)
%xlabel('Time in sec')
%ylabel('X body accel in g')
%subplot(212)
%plot(t0,y_a_B/g)
%xlabel('Time in sec')
%ylabel('Y body accel in g')
%
%figure (10)
%subplot(211)
%plot(t0,x_v_B)
%xlabel('Time in sec')
%ylabel('X body velocity in m/sec')
%subplot(212)
%plot(t0,y_v_B)
%xlabel('Time in sec')
%ylabel('Y body velocity in m/sec')
% Define inertial sensor parameters "accelerate mpu9150"

% accelerate (bias)
ba0=1*0.1*g;                             % initial accel bias in g in along-direction
bp0=-0.1*g;                              % initial accel bias in g in perpenticular-direction
err_factor = 1.0;
sig_xr_0 = err_factor*0.01*g/3600;         % accel bias stability in g/sec-sqrt(sec) in along-direction
sig_yr_0 = err_factor*0.02*g/3600;         % accel bias stability in g/sec-sqrt(sec) in penpenticular-direction
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
sig_x=opticalflow_err_factor*0.01;              % optical flow measurement noise in meter/sec
sig_y=opticalflow_err_factor*0.02;              % optical flow measurement noise in meter/sec in y-direction
%
sig_x=opticalflow_err_factor*0.1;              % optical flow measurement noise in meter/sec
sig_y=opticalflow_err_factor*0.2;              % optical flow measurement noise in meter/sec in y-direction
[vxm,vym,nvx,nvy]=XYvm(x_v_B,y_v_B,m,n,sig_x,sig_y);%含有bias optimal flow
% calculator "Q" use bias & noise
sig_bx=sig_bx_0;%(noise)
sig_by=sig_by_0;%(noise)
sig_xr=sig_xr_0;%(bias)
sig_yr=sig_yr_0;%(bias)

delta_t = dt;                           % delta time for simulating the true dynamics = 0.01 sec
%delta_s = delta_t;                      % sampling at every 0.01 second for the Kalman filter
delta_s = 10*delta_t;                      % sampling at every 0.01 second for the Kalman filter
[sensor_step,propagation_step]=propagate_step(T,delta_t,delta_s);

% Define the initial conditions for the inertial sensors
% Introduce +- 0.1m/sec initial errors in Vx_h and Vy_h; 0.1 deg error in
% psi_h
Vxerr = 0.1;
Vyerr = -0.1;
%===============================================================
% 2-state kalman filters
[axm_h,aym_h,Vx_h,Vy_h,vxm_h,vym_h,bx_h,by_h]=initial_estimate_value1(m,vxm,vym,axm,aym,x_v_B,y_v_B,Vxerr,Vyerr);

% Define the initial conditions for the Kalman Filter
[x_h,y_h,P00_x,P00_y]=define_initial_condition1(Vxerr,Vyerr,ba0,bp0);

% Po1(1)=sqrt(P(1,1)); Po2(1)=sqrt(P(2,2)); Po3(1)=sqrt(P(3,3));
% Po4(1)=sqrt(P(4,4)); Po5(1)=sqrt(P(5,5)); Po6(1)=sqrt(P(6,6));
% Po7(1)=sqrt(P(7,7)); Po8(1)=sqrt(P(8,8));
% Pn1(1)=sqrt(P(1,1)); Pn2(1)=sqrt(P(2,2)); Pn3(1)=sqrt(P(3,3));
% Pn4(1)=sqrt(P(4,4)); Pn5(1)=sqrt(P(5,5)); Pn6(1)=sqrt(P(6,6));
% Pn7(1)=sqrt(P(7,7)); Pn8(1)=sqrt(P(8,8));

% Define the constant matrices
H_x = zeros(1,2);H_x(1,1) = 1;
H_y = zeros(1,2);H_y(1,1) = 1;
R_x = sig_x^2;
R_y = sig_y^2;
F = [0 -1
     0  0]; 

% ============================================================
% Start the simulation run
k=2;
% ============================================================

for i=1:sensor_step
	for j=1:propagation_step
        k=1+j+(i-1)*propagation_step;
        
        bx_h(k)=bx_h(k-1);
        by_h(k)=by_h(k-1);
        %
        [axm_h,aym_h,Vx_h,Vy_h]=inertial_navigation_computation1(axm_h,aym_h,Vx_h,Vy_h,axm,aym,bx_h,by_h,k,dt);
        % 計算出phi和Q
        [phi,Q_x,Q_y]=define_Dymamic_equation1(F,sig_bx,sig_by,sig_xr,sig_yr,dt);
        % 卡爾曼的估計 
        [x_h,y_h,P00_x,P00_y]=Kalman_Filter_estimate1(x_h,y_h,phi,P00_x,P00_y,Q_x,Q_y,dt);
    end                 % end of filter propagation step    
% 
    [x_h,y_h,P00_x,P00_y,K_x,K_y,x_update,y_update]=Kalman_Filter_update1(x_h,y_h,P00_x,P00_y,Vx_h,Vy_h,vxm,vym,H_x,H_y,R_x,R_y,k);
%    s_x_update(:,i) = x_update;
%    s_x_update(:,i) = y_update;

%     P_new=P;
%     [Pn1,Pn2,Pn3,Pn4,Pn5,Pn6,Pn7,Pn8]=save_P(P_new,(i+1),Pn1,Pn2,Pn3,Pn4,Pn5,Pn6,Pn7,Pn8);

    [Vx_h,Vy_h,bx_h,by_h]=upon_optimalflow_measurement(x_h,y_h,k,Vx_h,Vy_h,bx_h,by_h);
%    old_x_h = x_h;
%    old_y_h = y_h;
      
    % reset the error
    x_h=zeros(2,1);
    y_h=zeros(2,1);
end                     % end of one Monte Caro run

end


% ============================================================
% Plot the performance
% ============================================================
e_point = 1000;
n1 = 400;
n2 = 100;
n2 = length(t0);

figure (3)
plot(t0,x_p_N,t0,y_p_N)
xlabel('Time in sec')
ylabel('Position in m')
%
figure (4)
subplot(211)
plot(t0,x_v_B,t0,Vx_h)
%xlabel('Time in sec')
ylabel('X Velocity in m/sec')
subplot(212)
plot(t0,y_v_B,t0,Vy_h)
xlabel('Time in sec')
ylabel('Y Velocity in m/sec')
%
figure (5)
subplot(211)
plot(t0,x_a_B/g,t0,axm_h/g)
%xlabel('Time in sec')
ylabel('X accel in g')
subplot(212)
plot(t0,y_a_B/g,t0,aym_h/g)
%xlabel('Time in sec')
ylabel('Y accel in g')
%
figure (8)
subplot(211)
plot(t0,x_a_N/g,t0,x_a_B/g)
%xlabel('Time in sec')
ylabel('X accel in g')
subplot(212)
plot(t0,y_a_N/g,t0,y_a_B/g)
%xlabel('Time in sec')
ylabel('Y accel in g')
%
figure (6)
subplot(211)
plot(t0(n1:n2),x_v_B(n1:n2)-Vx_h(n1:n2)')
%xlabel('Time in sec')
ylabel('X velocity est err in m/sec')
subplot(212)
plot(t0(n1:n2),y_v_B(n1:n2)-Vy_h(n1:n2)')
%xlabel('Time in sec')
ylabel('Y velocity est err in m/sec')
%
figure(7);
subplot(2,1,1), plot(t0,bx/g,'r',t0,bx_h/g,'g');         
xlabel('Time in seconds');ylabel('Accel bias along-axis in g');
subplot(2,1,2), plot(t0,by/g,'r',t0,by_h/g,'g');         
xlabel('Time in seconds');ylabel('Accel bias perp-axis in g');

%
figure (1)
subplot(2,1,1), plot(t0,bx/g,'r');         
xlabel('Time in seconds');ylabel('Accel bias along-axis in g');
subplot(2,1,2), plot(t0,by/g,'r');         
xlabel('Time in seconds');ylabel('Accel bias perp-axis in g');
%
figure (2)
subplot(2,1,1), plot(t0(n1:n2),(bx(n1:n2)-bx_h(n1:n2))/g,'r');         
xlabel('Time in seconds');ylabel('Accel bias err in x-axis in g');
grid
subplot(2,1,2), plot(t0(n1:n2),(by(n1:n2)-by_h(n1:n2))/g,'r');         
xlabel('Time in seconds');ylabel('Accel bias perp-axis in g');
grid