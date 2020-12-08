function [R1m,R2m,R3m,R4m,nvx_r,nvy_r] = radio_sensor_m_4(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,x_p_N,y_p_N,h_1,h_2,h_3,h_4,sig_x_r,sig_y_r,n,m)
%radiosensor_err_factor = 1.0; R1m(i) = sqrt((xr1-x_p_N(i))^2+(yr1-y_p_N(i))^2+h^2)+nvx_r(i);
%sig_x_r=radiosensor_err_factor*0.1;              % radio sensor measurement noise in meter/sec
%sig_y_r=radiosensor_err_factor*0.2;              % radio sensor measurement noise in meter/sec in y-direction
nvx_r=normrnd(0,sig_x_r,1,n);
nvy_r=normrnd(0,sig_y_r,1,n);
R1m=zeros(m);
R2m=zeros(m);
R3m=zeros(m);
R4m=zeros(m);
for i=1:n
    R1m(i) = sqrt((xr1-x_p_N(i))^2+(yr1-y_p_N(i))^2+h_1^2)+nvx_r(i);
    R2m(i) = sqrt((xr2-x_p_N(i))^2+(yr2-y_p_N(i))^2+h_2^2)+nvy_r(i);
    R3m(i) = sqrt((xr3-x_p_N(i))^2+(yr3-y_p_N(i))^2+h_3^2)+nvx_r(i);
    R4m(i) = sqrt((xr4-x_p_N(i))^2+(yr4-y_p_N(i))^2+h_4^2)+nvy_r(i);
end
end