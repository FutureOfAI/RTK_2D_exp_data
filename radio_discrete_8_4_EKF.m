function[H,R,R1m_h,R2m_h,R3m_h,R4m_h]=radio_discrete_8_4_EKF(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,h_1,h_2,h_3,h_4,k,format_1)
if (format_1 == 1) % straight line
R1m_h(k) = sqrt((xr1-xpm_Nh(k))^2+(yr1-ypm_Nh(k))^2+(h_1-1.32)^2); 
R2m_h(k) = sqrt((xr2-xpm_Nh(k))^2+(yr2-ypm_Nh(k))^2+(h_2-1.32)^2);
R3m_h(k) = sqrt((xr3-xpm_Nh(k))^2+(yr3-ypm_Nh(k))^2+(h_3-1.32)^2);
R4m_h(k) = sqrt((xr4-xpm_Nh(k))^2+(yr4-ypm_Nh(k))^2+(h_4-1.32)^2);

elseif (format_1 == 2)% square
R1m_h(k) = sqrt((xr1-xpm_Nh(k))^2+(yr1-ypm_Nh(k))^2+(h_1-1.2^2)); 
R2m_h(k) = sqrt((xr2-xpm_Nh(k))^2+(yr2-ypm_Nh(k))^2+(h_2-1.2)^2);
R3m_h(k) = sqrt((xr3-xpm_Nh(k))^2+(yr3-ypm_Nh(k))^2+(h_3-1.2)^2);
R4m_h(k) = sqrt((xr4-xpm_Nh(k))^2+(yr4-ypm_Nh(k))^2+(h_4-1.2)^2);
end
r1_partial_x =-(xr1-xpm_Nh(k))/R1m_h(k);
r1_partial_y =-(yr1-ypm_Nh(k))/R1m_h(k);
r2_partial_x =-(xr2-xpm_Nh(k))/R2m_h(k); 
r2_partial_y =-(yr2-ypm_Nh(k))/R2m_h(k);
r3_partial_x =-(xr3-xpm_Nh(k))/R3m_h(k);
r3_partial_y =-(yr3-ypm_Nh(k))/R3m_h(k);
r4_partial_x =-(xr4-xpm_Nh(k))/R4m_h(k); 
r4_partial_y =-(yr4-ypm_Nh(k))/R4m_h(k);

H = [r1_partial_x 0 0 r1_partial_y 0 0 0 0
              r2_partial_x 0 0 r2_partial_y 0 0 0 0
              r3_partial_x 0 0 r3_partial_y 0 0 0 0
              r4_partial_x 0 0 r4_partial_y 0 0 0 0 ];
% if k < 1500
% R = 1*[sig_x_r^2 0 0 0
%     0 sig_y_r^2 0 0
%     0 0 sig_x_r^2 0
%     0 0 0 sig_y_r^2];   
% else
R = 100*[sig_x_r^2 0 0 0
    0 sig_y_r^2 0 0
    0 0 sig_x_r^2 0
    0 0 0 sig_y_r^2];  
% end

end