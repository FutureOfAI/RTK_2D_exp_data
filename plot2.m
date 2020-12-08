n1 = 2000;
n2 = 100;
n2 = k-1;
%
figure (9)
subplot(211)
plot(t0(n1:n2),psi(n1:n2)*r2d,t0(n1:n2),psi_h(n1:n2)*r2d)
xlabel('Time in seconds');ylabel('Angle in deg');
subplot(212)
plot(t0(n1:n2),wzm_h(n1:n2)*r2d,t0(n1:n2),wz(n1:n2)*r2d)
ylabel('Angle rate in deg/sec') 
%
figure (11)
subplot(211)
plot(t0(n1:n2),x_p_N(n1:n2)-xpm_Nh(n1:n2)')
ylabel('X err position in meter')
subplot(212)
plot(t0(n1:n2),y_p_N(n1:n2)-ypm_Nh(n1:n2)')
ylabel('Y err position in meter')
xlabel('Time in seconds')
%
figure (12)
subplot(211)
plot(t0(n1:n2),(psi(n1:n2)-psi_h(n1:n2)')*r2d)
xlabel('Time in seconds');ylabel('Angle error in deg');
subplot(212)
plot(t0(n1:n2),(wzm_h(n1:n2)-wz(n1:n2)')*r2d)
ylabel('Angle rate error in deg/sec') 
%
figure (10)
subplot(211)
plot(t0(n1:n2),(psi(n1:n2)-psi_h(n1:n2)')*r2d)
xlabel('Time in seconds');ylabel('Angle error in deg');
subplot(212)
plot(t0(n1:n2),(bz(n1:n2)-bz_h(n1:n2))*r2d)
xlabel('Time in seconds');ylabel('Gyro bias error in deg/sec');

