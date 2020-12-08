% program plot5.m - plot the flight profile
figure (1)
plot(time,roll_angle*r2d)
xlabel('Time in seconds')
ylabel('Roll angle in deg')
grid
%
figure (2)
plot(time,yaw_rate*r2d)
xlabel('Time in seconds')
ylabel('yaw angle rate in deg/sec')
grid
%
figure (3)
plot(time,yaw_angle*r2d)
xlabel('Time in seconds')
ylabel('yaw angle in deg')
grid
%
figure (4)
subplot(211)
plot(time,x_v_N)
ylabel('Velocity in x-axis, m/sec')
grid
subplot(212)
plot(time,y_v_N)
xlabel('Time in seconds')
ylabel('Velocity in y-axis, m/sec')
grid
%
figure (5)
plot(x_p_N,y_p_N)
xlabel('Platform x-axis position in m')
ylabel('Platform y-axis position in m')
grid
%
figure (6)
subplot(211)
plot(time,roll_angle*r2d)
ylabel('Platform roll angle in deg')
grid
subplot(212)
plot(time,yaw_angle*r2d)
xlabel('Time in seconds')
ylabel('Platform yaw angle in deg')
grid
