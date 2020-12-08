figure (3)
subplot(311)
plot(x_p_N,y_p_N,coordinate(:,2),coordinate(:,3),'r*')
xlabel('X position in m')
ylabel('Y position in m')
grid
subplot(312)
plot(t00,x_v_N,'r',t00,y_v_N)
xlabel('Time in sec')
ylabel('Velocity in m/s')
grid
subplot(313)
plot(t0,x_a_B/g,'r',t0,y_a_B/g,'g')
ylabel('X-Y body accel in g')
xlabel('Time in sec')
grid