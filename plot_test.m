n1 = 400;
n2 = k-1;

figure (17)
% plot(t0(n1:n2),(psi(n1:n2)-psi_h(n1:n2)')*r2d)
plot(t0(n1:n2),(psi_h(n1:n2)')*r2d)
xlabel('Time in seconds');ylabel('Psi in deg');
ylabel('Angle in deg') 