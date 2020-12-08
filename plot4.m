% plot4.m
%
f = [0 10 20 50 100];           % Kalman filter update rate in Hz
%
fx6c2 = [0.3 0.16047 0.12251 0.087141 0.06726]/3*100;
fy6c2 = [0.3 0.15619 0.11887 0.084154 0.064643]/3*100;
%
fx6s2 = [0.3 0.27 0.2446 0.1782 0.1436]/3*100;
fy6s2 = [0.3 0.1350 0.109 0.0710 0.0542]/3*100;
%
fx6c4 = [0.3/sqrt(2) 0.14609 0.10474 0.074619 0.057845]/3*100;
fy6c4 = [0.3/sqrt(2) 0.14691 0.10693 0.079831 0.062765]/3*100;
%
fx6s4 = [0.3/sqrt(2) 0.12512 0.085906 0.064152 0.049564]/3*100;
fy6s4 = [0.3/sqrt(2) 0.12508 0.085642 0.064295 0.04972]/3*100;
%
figure (1)
plot(f,fx6c2,'r',f,fy6c2,'b',f,fx6c4,'k-',f,fy6c4,'g-')
xlabel('Filter update rate in Hz')
ylabel('One sigma position err in cm')
grid
axis([0 100 0 10])
%
figure (2)
plot(f,fx6s2,'r',f,fy6s2,'b',f,fx6s4,'k-',f,fy6s4,'g-')
xlabel('Filter update rate in Hz')
ylabel('One sigma position err in cm')
grid
axis([0 100 0 10])
%
figure (3)
plot(f,fx6s2,'r',f,fy6s2,'b',f,fx6c2,'k-',f,fy6c2,'g-')
xlabel('Filter update rate in Hz')
ylabel('One sigma position err in cm')
grid
axis([0 100 0 10])
%
figure (4)
plot(f,fx6c2,'r',f,fy6c2,'b')
xlabel('Filter update rate in Hz')
ylabel('One sigma position err in cm')
grid
axis([0 100 0 10])