
n1 = 1;
n2 = 100;
n2 = k-1;
%t00 = time;
% for i = 1:n2
% xpm_Nh_1(i,:) = xpm_Nh(i,:);
% ypm_Nh_1(i,:) = ypm_Nh(i,:);
% end
for i = 1:count-9 %20s
mean_3(1,i) = true_trajectory-(xpm_Nh(i)-0.02);
mean_4(1,i) = (ypm_Nh(i)+0.01)-position(1,i);    %% y error 
end
%% mean , std , 3sigma  x position
% mean_pis_h = mean(psi_h(1:count-9)) * 180/pi;
mean_KF = mean(mean_3(1,2000:count-9));
std_KF = std(mean_3(1,2000:count-9));
% % bias_err+3*std
% sigma_KF = abs(mean_KF) + (3*std_KF);
% %% mean , std , 3sigma  y position
% mean_LS_y = mean(mean_4(1,2000:count-9));
% std_LS_y = std(mean_4(1,2000:count-9));
% % bias_err+3*std
% sigma_LS_y = abs(mean_LS_y) + (3*std_LS_y);

% mean_Mu_z_1 = mean(Mu_z_1(1000:1995,1));
% mean_Mu_z_2 = mean(Mu_z_2(1000:1995,1));
% mean_Mu_z_3 = mean(Mu_z_3(1000:1995,1));
% mean_Mu_z_4 = mean(Mu_z_4(1000:1995,1));
% std_Mu_z_1 = std(Mu_z_1(1000:1995,1));
% std_Mu_z_2 = std(Mu_z_2(1000:1995,1));
% std_Mu_z_3 = std(Mu_z_3(1000:1995,1));
% std_Mu_z_4 = std(Mu_z_4(1000:1995,1));
% sigma_Mu_z_1 = abs(mean_Mu_z_1) + (3*std_Mu_z_1);
% sigma_Mu_z_2 = abs(mean_Mu_z_2) + (3*std_Mu_z_2);
% sigma_Mu_z_3 = abs(mean_Mu_z_3) + (3*std_Mu_z_3);
% sigma_Mu_z_4 = abs(mean_Mu_z_4) + (3*std_Mu_z_4);

%% CDF (Cumulative distribution functions)


% pd_KF = makedist('Normal',mean_KF,std_KF);%('Normal',mean,std)
% cdf_KF = cdf(pd_KF,ranger);


figure (1)
plot(y,x,'--r',position(1,:),position(2,:),y_anchor_1,x_anchor_1,'r*')
hold on
% plot(5.4,1.4,'diamond',5.4,9.2,'rsquare',5.4,1.4,'bx','MarkerSize',13,'LineWidth',2) % 10 ¨Ó¦^
% plot(5.4,0.74,'diamond',5.4,6.8,'rsquare','MarkerSize',13,'LineWidth',2)%,5.4,0.74,'bx'
plot(0.85,3.35,'diamond',5.45,4.11,'rsquare','MarkerSize',13,'LineWidth',2)%,5.4,0.74,'bx'
legend('Reference path','Least square','Anchor','Start point','End point');%'Reentry point',
hold on
plot(y_anchor_2,x_anchor_2,'r*',y_anchor_3,x_anchor_3,'r*',y_anchor_4,x_anchor_4,'r*')
title('Least square')
xlabel('Y position in m')
ylabel('X position in m')
axis([0 10 0 10])
set(gca,'YDir','reverse')
grid

figure (2)
plot(y,x,'--r',position_2(1,:),position_2(2,:),y_anchor_1,x_anchor_1,'r*')
hold on
plot(5.4,0.74,'diamond',5.4,6.8,'rsquare',5.4,0.74,'bx','MarkerSize',13,'LineWidth',2)
legend('Reference path','Least square','Anchor','Start point','Reentry point','End point');
hold on
plot(y_anchor_2,x_anchor_2,'r*',y_anchor_3,x_anchor_3,'r*',y_anchor_4,x_anchor_4,'r*')
anchor_d=p;
title(['Trilateration method not Anchor',num2str(anchor_d)])
xlabel('Y position in m')
ylabel('X position in m')
axis([0 10 0 10])
%axis([2.5 5.5 0.5 4])
set(gca,'YDir','reverse')
grid
    
%% plot line
line([1,10],[1.12 1.12])

% y_1 = [3.225 3.225];
% x_1 = [1 4];    
% 
% y_2 = [5.153 5.153];
% x_2 = [1 4];    

% y_3 = [3  5.5];
% x_3 = [1.5 1.5]; 
% 
% y_4 = [3 5.5];
% x_4 = [3.395 3.395]; 

y_1 = [3.225 3.225];
x_1 = [1.5 3.395];    

y_2 = [5.153 5.153];
x_2 = [1.5 3.395];

y_3 = [3.225  5.153];
x_3 = [1.5 1.5]; 

y_4 = [3.225 5.153];
x_4 = [3.395 3.395]; 

y_5 = [6.3 7.1];
x_5 = [4.64 4.64]; 

figure (3)

%% paper test 1 line straight 
% plot(position_2(1,:),position_2(2,:),'g.',ypm_Nh(1:3991),xpm_Nh(1:3991),'r')%
% hold on
% legend('TRI','EKF');
% hold on
%  
% plot(ypm_Nh(1),xpm_Nh(1),'diamond',ypm_Nh(3991),xpm_Nh(3991),'bx','MarkerSize',13,'LineWidth',2)

% % paper use test
% plot(position(1,:),position(2,:),'MarkerSize',13,'LineWidth',1)%position_2(1,:),position_2(2,:),'g.',ypm_Nh(1:3991),xpm_Nh(1:3991)
% hold on
% 
%  plot(position(1,1),position(2,1),'rdiamond',position(1,3991),position(2,3991),'rx',y_5,x_5,'--r','MarkerSize',13,'LineWidth',2)
% % plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(3991),xpm_Nh(3991),'rx','MarkerSize',13,'LineWidth',2)
% 
% 
% legend('LS','Start point','End point','Reference path');
% % legend('EKF','Start point','End point');


%ypm_Nh(2000)+0.01,xpm_Nh(2000)-0.02,'rsquare',
% legend('Start point','Reentry point','End point','Reference path');%'Reference path','Least square','Anchor',
% plot(ypm_Nh_1(1:2000),xpm_Nh_1(1:2000),'g',ypm_Nh(1:2000),xpm_Nh(1:2000),'b')%ypm_Nh_1(1:2000),xpm_Nh_1(1:2000)
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(3991),xpm_Nh(3991),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% hold on
% plot(y_3,x_3,'--r')
% hold on
% plot(yr1,xr1,'r*')
% hold on
% legend('Start point','End point');%'Reference path','Least square','Anchor',
%% paper square ekf_ls
% plot(position(1,1:2000),position(2,1:2000),'g',ypm_Nh(1:2000),xpm_Nh(1:2000),'b')
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(2000),xpm_Nh(2000),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% plot(y_3,x_3,'--r')
% hold on
% legend('LS','EKF','Start point','End point','Reference path');%'Reference path','Least square','Anchor',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% hold on
% plot(y_1,x_1,'--r',y_2,x_2,'--r',y_4,x_4,'--r',y_5,x_5,'--r')
% hold on

%%   (5/17) Square
% plot(ypm_Nh(1:1900),xpm_Nh(1:1900),'MarkerSize',9,'LineWidth',1)%position_2(1,1:2000),position_2(2,1:2000)
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(1900),xpm_Nh(1900),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% plot(y_3,x_3,'--r','MarkerSize',9,'LineWidth',1)
% hold on
% legend('TRI','Start point','End point','Reference path');%'Reference path','Least square','Anchor','LS',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*','MarkerSize',9,'LineWidth',1)
% hold on
% plot(y_1,x_1,'--r',y_2,x_2,'--r',y_4,x_4,'--r',y_5,x_5,'--r','MarkerSize',9,'LineWidth',1)
% hold on

%% 
% plot(ypm_Nh(1:3991),xpm_Nh(1:3991),'MarkerSize',9,'LineWidth',1)%position_2(1,1:2000),position_2(2,1:2000)
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(3991),xpm_Nh(3991),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% plot(y_3,x_3,'--r','MarkerSize',9,'LineWidth',1)
% hold on
% legend('TRI','Start point','End point','Reference path');%'Reference path','Least square','Anchor','LS',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*','MarkerSize',9,'LineWidth',1)
% hold on
% plot(y_1,x_1,'--r',y_2,x_2,'--r',y_4,x_4,'--r',y_5,x_5,'--r','MarkerSize',9,'LineWidth',1)
% hold on



%%   (5/17) circle
% plot(position_2(1,1:1350),position_2(2,1:1350),'MarkerSize',9,'LineWidth',1)%ypm_Nh(1:1350),xpm_Nh(1:1350)
% hold on
% plot(position_2(1,1),position_2(2,1),'rdiamond',position_2(1,1350),position_2(2,1350),'rx','MarkerSize',9,'LineWidth',2)%position_2(1,1),position_2(2,1)
% hold on
% plot(x_circle_1,y_circle_1,'r--','MarkerSize',9,'LineWidth',1)
% hold on
% legend('TRI','Start point','End point','Reference path');%'Reference path','Least square','Anchor','LS',
% hold on



%% square c.05_c2
% plot(ypm_Nh(1:2000),xpm_Nh(1:2000),'b',ypm_Nh_1(1:2000),xpm_Nh_1(1:2000),'g','MarkerSize',9,'LineWidth',1)
% hold on
% plot(y_3,x_3,'--r','MarkerSize',9,'LineWidth',1)
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(2000),xpm_Nh(2000),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% legend('c = 2','c = 0.5','Reference path','Start point','End point');%'Reference path','Least square','Anchor',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% hold on
% plot(y_1,x_1,'--r',y_2,x_2,'--r',y_4,x_4,'--r',y_5,x_5,'--r','MarkerSize',9,'LineWidth',1)
% hold on
%% other trajectory
% plot(position(1,:),position(2,:),'g',ypm_Nh(1:3991),xpm_Nh(1:3991),'b')
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(3991),xpm_Nh(3991),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% plot(yr1,xr1,'r*')
% hold on
% legend('LS','EKF','Start point','End point','Anchor');%'Reference path','Least square','Anchor',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% hold on

% %% other trajectory 5/17
% plot(ypm_Nh(1:3991),xpm_Nh(1:3991),'MarkerSize',9,'LineWidth',1)%,ypm_Nh(1:3991),xpm_Nh(1:3991)'g',position_2(1,:),position_2(2,:)
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(3991),xpm_Nh(3991),'rx','MarkerSize',9,'LineWidth',2)%position_2(1,1),position_2(2,1)
% hold on
% plot(yr1,xr1,'r*','MarkerSize',12,'LineWidth',1)
% hold on
% legend('TRI','Start point','End point','Anchor');%'Reference path','Least square','Anchor',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*','MarkerSize',12,'LineWidth',1)
% hold on

%% 5/18
% plot(ypm_Nh(1:3991),xpm_Nh(1:3991),'b',ypm_Nh_1(1:3990),xpm_Nh_1(1:3990),'g','MarkerSize',9,'LineWidth',1)
% hold on
% 
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(3991),xpm_Nh(3991),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% plot(yr1,xr1,'r*')
% hold on
% legend('c = 2','c = 0.3','Start point','End point','Anchor');%'Reference path','Least square','Anchor',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% hold on
%%  point 5/17
% plot(position(1,:),position(2,:),'MarkerSize',9,'LineWidth',1)%,ypm_Nh(1:3991),xpm_Nh(1:3991)'g',position_2(1,:),position_2(2,:)
% hold on
% plot(position(1,1),position(2,1),'rdiamond',position(1,3991),position(2,3991),'rx','MarkerSize',9,'LineWidth',2)%position_2(1,1),position_2(2,1)
% hold on
% legend('LS','Start point','End point','Anchor');%'Reference path','Least square','Anchor',
% hold on


%% circle EKF
% plot(ypm_Nh(1:1350),xpm_Nh(1:1350),'b',ypm_Nh_1(1:1350),xpm_Nh_1(1:1350),'g',x_circle_1,y_circle_1,'r--','MarkerSize',9,'LineWidth',1)
% hold on
% plot(ypm_Nh(1),xpm_Nh(1),'rdiamond',ypm_Nh(1350),xpm_Nh(1350),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% legend('c = 2','c = 0.1','Reference path','Start point','End point');%'Reference path','Least square','Anchor',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% hold on

%% circle TRI & LS
% plot(position_2(1,1:1350),position_2(2,1:1350),'g',x_circle_1,y_circle_1,'r--')
% hold on
% plot(position_2(1,1),position_2(2,1),'rdiamond',position_2(1,1350),position_2(2,1350),'rx','MarkerSize',9,'LineWidth',2)
% hold on
% legend('TRI','Reference path','Start point','End point');%'Reference path','Least square','Anchor',
% hold on
% plot(yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% hold on


set(gca,'YDir','reverse') 
xlabel('Y position in m')
ylabel('X position in m')
%  axis([0 12 0 10])
 axis([6 7.5 4.3 5])
% axis([2.5 6 0.5 4.5])
% axis([6.65 6.85 4.50 4.8])
grid


hold on

%% other trajectory
% % figure (3)
% % plot(ypm_Nh(1),xpm_Nh(1),'diamond',ypm_Nh(2000),xpm_Nh(2000),'rsquare',ypm_Nh(3991),xpm_Nh(3991),'bx','MarkerSize',13,'LineWidth',2)
% % hold on
% % % plot(y_3,x_3,'--r')
% % % legend('Start point','Reentry point','End point','Reference path');%'Reference path','Least square','Anchor',
% % legend('Start point','Reentry point','End point');%'Reference path','Least square','Anchor',
% % hold on
% % % plot(position(1,1:500),position(2,1:500),'b',position(1,500:1000),position(2,500:1000),'r',position(1,1000:1500),position(2,1000:1500),'g',position(1,1500:2000),position(2,1500:2000),'m',position(1,2000:2500),position(2,2000:2500),'c',position(1,2500:3000),position(2,2500:3000),'k',position(1,3000:3500),position(2,3000:3500),'y',position(1,3500:3991),position(2,3500:3991),'k',yr1,xr1,'r*',yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% % plot(ypm_Nh(1:500),xpm_Nh(1:500),'b',ypm_Nh(500:1000),xpm_Nh(500:1000),'r',ypm_Nh(1000:1500),xpm_Nh(1000:1500),'g',ypm_Nh(1500:2000),xpm_Nh(1500:2000),'m',ypm_Nh(2000:2500),xpm_Nh(2000:2500),'c',ypm_Nh(2500:3000),xpm_Nh(2500:3000),'k',ypm_Nh(3000:3500),xpm_Nh(3000:3500),'y',ypm_Nh(3500:3991),xpm_Nh(3500:3991),'k',yr1,xr1,'r*',yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% % hold on
% % % plot(y_1,x_1,'--r',y_2,x_2,'--r',y_4,x_4,'--r')
% % % hold on
% % % plot(position_2(1,n1:n2),position_2(2,n1:n2),'g',ypm_Nh(n1:n2),xpm_Nh(n1:n2),'r',yr1,xr1,'r*',yr2,xr2,'r*',yr3,xr3,'r*',yr4,xr4,'r*')
% % % ,ypm_Nh(n1:n2),xpm_Nh(n1:n2),'r',position_2(1,:),position_2(2,:),'g'
% % set(gca,'YDir','reverse') 
% % xlabel('Y position in m')
% % ylabel('X position in m')
% % axis([0 12 0 10])
% % % axis([6 7.5 4.3 5])
% % % axis([2.5 6 0.5 4.5])
% % grid
% % hold on

figure (17)
% subplot(211)
% % plot(t0(n1:n2),((low_gro(3,n1:n2)'-bz_h(n1:n2))*r2d))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
% plot(t0(n1:n2),((gro(3,n1:n2)'-bz_h(n1:n2))*r2d))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
% xlabel('Time in seconds');ylabel('Angular rate error in deg/sec');
% subplot(212)
% plot(t0(n1:n2),((axm(n1:n2)/g)),t0(n1:n2),((low_acc(2,n1:n2)')))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
plot(t0(n1:n2),((acc(2,n1:n2)')))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time(s)'); ylabel('aym in g');%t0(n1:n2),((axm(n1:n2)/g)),

figure (19)
% subplot(311)
plot(t0(n1:n2),(((acc(1,n1:n2)'))))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
% plot(t0(n1:n2),((aym(n1:n2)/g)),t0(n1:n2),((low_acc(1,n1:n2)')))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time(s)'); ylabel('axm in g');%(aym(n1:n2)/g)),t0(n1:n2),

figure (20)
subplot(411)
plot(t0(n1:n2),pos(1,n1:n2),'b')%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time in seconds'); ylabel('distance in m');
subplot(412)
plot(t0(n1:n2),pos(2,n1:n2),'b')%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time in seconds'); ylabel('distance in m');
subplot(413)
plot(t0(n1:n2),pos(3,n1:n2),'b')%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time in seconds'); ylabel('distance in m');
subplot(414)
plot(t0(n1:n2),pos(4,n1:n2),'b')%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time in seconds'); ylabel('distance in m');

figure(9);
subplot(3,1,1), plot(t0,bx/g,'r',t0,bx_h,'g');         
xlabel('Time in seconds');ylabel('Accel bias along-axis in g');
subplot(3,1,2), plot(t0,by/g,'r',t0,by_h,'g');         
xlabel('Time in seconds');ylabel('Accel bias perp-axis in g');
subplot(3,1,3), plot(t0,bz*r2d,'r',t0,bz_h*r2d,'g');
xlabel('Time in seconds');ylabel('Gyro bias in deg');
%
figure (10)
subplot(3,1,1), plot(t0(n1:n2),bx_h(n1:n2),'b');
hold on
plot(t0(1263),bx_h(1263),'r+','MarkerSize',15);
xlabel('Time in seconds');ylabel('Accel bias along x-axis in g');
subplot(3,1,2), plot(t0(n1:n2),by_h(n1:n2),'b');   
hold on
plot(t0(1263),by_h(1263),'r+','MarkerSize',15);
xlabel('Time in seconds');ylabel('Accel bias along y-axis in g');
subplot(3,1,3), 
plot(t0(n1:n2),bz_h(n1:n2)*r2d,'b');
hold on
plot(t0(1263),bz_h(1263)*r2d,'r+','MarkerSize',15);
xlabel('Time in seconds');ylabel('Gyro bias in deg/sec');
hold off
%
figure (11)
subplot(2,1,1), plot(t0(n1:n2),(bx(n1:n2)-bx_h(n1:n2)),'r');         
xlabel('Time in seconds');ylabel('Accel bias err in x-axis in g');
grid
subplot(2,1,2), plot(t0(n1:n2),(by(n1:n2)-by_h(n1:n2)),'r');         
xlabel('Time in seconds');ylabel('Accel bias y-axis in g');
grid

figure (12)
subplot(211)
plot(t0(n1:n2),xam_Nh(n1:n2))
xlabel('Time in seconds')
ylabel('X accelerometer in g')
subplot(212)
plot(t0(n1:n2),yam_Nh(n1:n2))
xlabel('Time in seconds')
ylabel('Y accelerometer in g')

figure (13)
% subplot(311)
subplot(211)
plot(t0(n1:n2),xvm_Nh(n1:n2))
xlabel('Time in seconds')
ylabel('X velocity in m/sec')
subplot(212)
plot(t0(n1:n2),yvm_Nh(n1:n2))
xlabel('Time in seconds')
ylabel('Y velocity in m/sec')
%% line straight
% figure (14)
% % subplot(311)
% subplot(211)
% plot(t0(n1:n2),position_2(2,n1:n2),'MarkerSize',13,'LineWidth',1)%xpm_Nh(n1:n2)
% xlabel('Time in seconds')
% ylabel('X positioning in m')
% axis([0 40 4.5 4.9])
% legend('X');
% subplot(212)
% plot(t0(n1:n2),position_2(1,n1:n2),'MarkerSize',13,'LineWidth',1)%ypm_Nh(n1:n2)
% xlabel('Time in seconds')
% ylabel('Y positioning in m')
% axis([0 40 6 7.5])
% legend('Y');

%% square
figure (14)
% subplot(311)
subplot(211)
plot(t0(n1:n2),position(2,n1:n2),'MarkerSize',13,'LineWidth',1)%position_2(2,n1:n2)
xlabel('Time in seconds')
ylabel('X positioning in m')
% axis([0 13.5 -1 6])
% axis([0 40 0 10])
% axis([0 40 4.53 4.8])

legend('X');
subplot(212)
plot(t0(n1:n2),position(1,n1:n2),'MarkerSize',13,'LineWidth',1)%position_2(1,n1:n2)
xlabel('Time in seconds')
ylabel('Y positioning in m')
% axis([0 13.5 2 6])
% axis([0 40 3 7])
% axis([0 40 6.65 6.85])

legend('Y');

% figure(16)
% title('Cumulative distribution functions')
% plot(ranger,cdf_Tri,'--',ranger,cdf_LS,'-.',ranger,cdf_KF)%,1:1500,y__1
% legend('Tri','LS','KF')
% xlabel('position error[m]')
% ylabel('CDF')
% grid
%   Mu
% figure (15)
% % subplot(311)
% subplot(411)
% plot(t0(n1:2:n2),Mu_z_1(n1:n2/2),'b')
% axis([0 n2/100 -2 2])
% ylabel('\mu 1')
% xlabel('Time in seconds')
% subplot(412)
% plot(t0(n1:2:n2),Mu_z_2(n1:n2/2),'b')
% axis([0 n2/100 -2 2])
% ylabel('\mu 2')
% xlabel('Time in seconds')
% subplot(413)
% plot(t0(n1:2:n2),Mu_z_3(n1:n2/2),'b')
% axis([0 n2/100 -2 2])
% ylabel('\mu 3')
% xlabel('Time in seconds')
% subplot(414)
% plot(t0(n1:2:n2),Mu_z_4(n1:n2/2),'b')
% axis([0 n2/100 -2 2])
% ylabel('\mu 4')
% xlabel('Time in seconds')
% color
% figure (15)
% % subplot(311)
% subplot(411)
% plot(t0(1:2:500),Mu_z_1(1:500/2),'b',t0(500:2:1000),Mu_z_1(250:1000/2),'r',t0(1000:2:1500),Mu_z_1(500:1500/2),'g',t0(1500:2:2000),Mu_z_1(750:2000/2),'m',t0(2000:2:2500),Mu_z_1(1000:2500/2),'c',t0(2500:2:3000),Mu_z_1(1250:3000/2),'k',t0(3000:2:3500),Mu_z_1(1500:3500/2),'y',t0(3500:2:3990),Mu_z_1(1750:3990/2),'k')%
% axis([0 n2/100 -2 2])
% ylabel('\mu 1')
% xlabel('Time in seconds')
% subplot(412)
% plot(t0(1:2:500),Mu_z_2(1:500/2),'b',t0(500:2:1000),Mu_z_2(250:1000/2),'r',t0(1000:2:1500),Mu_z_2(500:1500/2),'g',t0(1500:2:2000),Mu_z_2(750:2000/2),'m',t0(2000:2:2500),Mu_z_2(1000:2500/2),'c',t0(2500:2:3000),Mu_z_2(1250:3000/2),'k',t0(3000:2:3500),Mu_z_2(1500:3500/2),'y',t0(3500:2:3990),Mu_z_2(1750:3990/2),'k')%
% axis([0 n2/100 -2 2])
% ylabel('\mu 2')
% xlabel('Time in seconds')
% subplot(413)
% plot(t0(1:2:500),Mu_z_3(1:500/2),'b',t0(500:2:1000),Mu_z_3(250:1000/2),'r',t0(1000:2:1500),Mu_z_3(500:1500/2),'g',t0(1500:2:2000),Mu_z_3(750:2000/2),'m',t0(2000:2:2500),Mu_z_3(1000:2500/2),'c',t0(2500:2:3000),Mu_z_3(1250:3000/2),'k',t0(3000:2:3500),Mu_z_3(1500:3500/2),'y',t0(3500:2:3990),Mu_z_3(1750:3990/2),'k')%
% axis([0 n2/100 -2 2])
% ylabel('\mu 3')
% xlabel('Time in seconds')
% subplot(414)
% plot(t0(1:2:500),Mu_z_4(1:500/2),'b',t0(500:2:1000),Mu_z_4(250:1000/2),'r',t0(1000:2:1500),Mu_z_4(500:1500/2),'g',t0(1500:2:2000),Mu_z_4(750:2000/2),'m',t0(2000:2:2500),Mu_z_4(1000:2500/2),'c',t0(2500:2:3000),Mu_z_4(1250:3000/2),'k',t0(3000:2:3500),Mu_z_4(1500:3500/2),'y',t0(3500:2:3990),Mu_z_4(1750:3990/2),'k')%
% axis([0 n2/100 -2 2])
% ylabel('\mu 4')
% xlabel('Time in seconds')
% %% not one anchor 


% figure (16)
% % subplot(311)
% subplot(511)
% plot(t0(n1:2:n2),Mu_count__1(1,n1:n2/2),'b')
% axis([0 n2/100 0 2])
% ylabel('Not Anchor 1')
% xlabel('Time in seconds')
% subplot(512)
% plot(t0(n1:2:n2),Mu_count__1(2,n1:n2/2),'b')
% axis([0 n2/100 0 2])
% ylabel('Not Anchor 2')
% xlabel('Time in seconds')
% subplot(513)
% plot(t0(n1:2:n2),Mu_count__1(3,n1:n2/2),'b')
% axis([0 n2/100 0 2])
% ylabel('Not Anchor 3')
% xlabel('Time in seconds')
% subplot(514)
% plot(t0(n1:2:n2),Mu_count__1(4,n1:n2/2),'b')
% axis([0 n2/100 0 2])
% ylabel('Not Anchor 4')
% xlabel('Time in seconds')
% subplot(515)
% plot(t0(n1:2:n2),Mu_count__1(5,n1:n2/2),'b')
% axis([0 n2/100 0 2])
% ylabel('All Anchor')
% xlabel('Time in seconds')
% % %% not two anchor 
% figure (17)
% % subplot(311)
% subplot(611)
% plot(t0(n1:2:n2),Mu_count__2(1,n1:n2/2),'b')% not anchor1 anchor2
% axis([0 n2/100 0 2])
% ylabel('Not Anchor1 & 2')
% xlabel('Time in seconds')
% subplot(612)
% plot(t0(n1:2:n2),Mu_count__2(2,n1:n2/2),'b')% not anchor1 anchor3
% axis([0 n2/100 0 2])
% ylabel('Not Anchor1 & 3')
% xlabel('Time in seconds')
% subplot(613)
% plot(t0(n1:2:n2),Mu_count__2(3,n1:n2/2),'b')% not anchor1 anchor4
% axis([0 n2/100 0 2])
% ylabel('Not Anchor1 & 4')
% xlabel('Time in seconds')
% subplot(614)
% plot(t0(n1:2:n2),Mu_count__2(4,n1:n2/2),'b')% not anchor2 anchor3
% axis([0 n2/100 0 2])
% ylabel('Not Anchor2 & 3')
% xlabel('Time in seconds')
% subplot(615)
% plot(t0(n1:2:n2),Mu_count__2(5,n1:n2/2),'b')% not anchor2 anchor4
% axis([0 n2/100 0 2])
% ylabel('Not Anchor2 & 4')
% xlabel('Time in seconds')
% subplot(616)
% plot(t0(n1:2:n2),Mu_count__2(6,n1:n2/2),'b')% not anchor3 anchor4
% axis([0 n2/100 0 2])
% ylabel('Not Anchor3 & 4')
% xlabel('Time in seconds')
% figure (18)
% % subplot(311)
% subplot(411)
% plot(t0(n1:2:n2),Mu_count__3(1,n1:n2/2),'b')% not anchor1 anchor2
% axis([0 n2/100 0 2])
% ylabel('Not Anchor1 & 2 & 3')
% xlabel('Time in seconds')
% subplot(412)
% plot(t0(n1:2:n2),Mu_count__3(2,n1:n2/2),'b')% not anchor1 anchor3
% axis([0 n2/100 0 2])
% ylabel('Not Anchor1 & 2 & 4')
% xlabel('Time in seconds')
% subplot(413)
% plot(t0(n1:2:n2),Mu_count__3(3,n1:n2/2),'b')% not anchor1 anchor4
% axis([0 n2/100 0 2])
% ylabel('Not Anchor1 & 3 & 4')
% xlabel('Time in seconds')
% subplot(414)
% plot(t0(n1:2:n2),Mu_count__3(4,n1:n2/2),'b')% not anchor2 anchor3
% axis([0 n2/100 0 2])
% ylabel('Not Anchor2 & 3 & 4')
% xlabel('Time in seconds')


