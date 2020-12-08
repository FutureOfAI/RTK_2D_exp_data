n1 = 1;
n2 = 500;
n2 = k-1;
%%  select 1 = original data ; select 2 = low pass filter data
select = 1;

if select == 1

figure (6)
plot(xpm_Nh(100:3990),ypm_Nh(100:3990),'b',coordinate(:,2),coordinate(:,3),'r*');
hold on
plot(xpm_Nh(1263),ypm_Nh(1263),'r+',xpm_Nh(2745),ypm_Nh(2745),'r+','MarkerSize',15);
text(coordinate(1,2),coordinate(1,3),' An1');
text(coordinate(2,2),coordinate(2,3),' An2');
text(coordinate(3,2),coordinate(3,3),' An3');
text(coordinate(4,2),coordinate(4,3),' An4');
text(coordinate(5,2),coordinate(5,3),' An5');
text(coordinate(6,2),coordinate(6,3),' An6');
text(coordinate(7,2),coordinate(7,3),' An7');
text(coordinate(8,2),coordinate(8,3),' An8');
text(coordinate(9,2),coordinate(9,3),' An9');
text(coordinate(10,2),coordinate(10,3),' An10');
% plot(x_p_N,y_p_N,'g*',xpm_Nh,ypm_Nh,'b',xr1,yr1,'r*',xr2,yr2,'r*',xr3,yr3,'r*',xr4,yr4,'r*')
xlabel('X position in m')
ylabel('Y position in m')
hold off
% figure (1)
% plot(t0(n1:n2),gro(3,n1:n2)*r2d)%
% xlabel('Time in seconds');ylabel('Angular rate in deg/sec');

% figure (2)
% plot(t0(n1:n2),((gro(3,n1:n2)'-bz_h(n1:n2))*r2d))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
% xlabel('Time in seconds');ylabel('Angular rate error in deg/sec');

% figure (3)
% plot(t0(n1:n2),(acc(1,n1:n2)'-bx_h(n1:n2)))
% xlabel('Time in seconds');ylabel('x accelerate error in g');

% figure (4)
% plot(t0(n1:n2),(acc(2,n1:n2)'-by_h(n1:n2)))
% xlabel('Time in seconds');ylabel('y accelerate error in g');

%square turn one circle is 1900 seconed
figure (5)
plot(t0(n1:n2),(psi_h(n1:n2)')*r2d)
xlabel('Time in seconds');ylabel('psi_h Angle in deg');

% figure (6)
% subplot(311)
% plot(t0(n1:n2),gro(3,n1:n2)*r2d)%
% xlabel('Time in seconds');ylabel('Angular rate in deg/sec');
% subplot(312)
% plot(t0(n1:n2),((gro(3,n1:n2)'-bz_h(n1:n2))*r2d))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
% xlabel('Time in seconds');ylabel('Angular rate error in deg/sec');
% subplot(313)
% plot(t0(n1:n2),(psi_h(n1:n2)')*r2d)
% xlabel('Time in seconds');ylabel('psi_h Angle in deg');

end
%% low pass filter
if select == 2

figure (1)
plot(t0(n1:n2),low_gro(3,n1:n2)*r2d)%
xlabel('Time in seconds');ylabel('Angular rate in deg/sec');

figure (2)
plot(t0(n1:n2),((low_gro(3,n1:n2)'-bz_h(n1:n2))*r2d))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time in seconds');ylabel('Angular rate error in deg/sec');

figure (3)
plot(t0(n1:n2),(acc(1,n1:n2)'-bx_h(n1:n2)))
xlabel('Time in seconds');ylabel('x accelerate error in g');

figure (4)
plot(t0(n1:n2),(acc(2,n1:n2)'-by_h(n1:n2)))
xlabel('Time in seconds');ylabel('y accelerate error in g');

figure (5)
plot(t0(n1:n2),(psi_h(n1:n2)')*r2d)
xlabel('Time in seconds');ylabel('psi_h Angle in deg');

figure (6)
subplot(311)
plot(t0(n1:n2),low_gro(3,n1:n2)*r2d)%
xlabel('Time in seconds');ylabel('Angular rate in deg/sec');
subplot(312)
plot(t0(n1:n2),((low_gro(3,n1:n2)'-bz_h(n1:n2))*r2d))%gro(3,n1:n2)'-bz_h(n1:n2))   linVelLP(n1:n2,3)
xlabel('Time in seconds');ylabel('Angular rate error in deg/sec');
subplot(313)
plot(t0(n1:n2),(psi_h(n1:n2)')*r2d)
xlabel('Time in seconds');ylabel('psi_h Angle in deg');

end