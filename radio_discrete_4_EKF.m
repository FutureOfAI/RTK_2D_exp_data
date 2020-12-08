function[H,R,R1m_h,R2m_h,R3m_h,R4m_h]=radio_discrete_4_EKF(xr1,yr1,xr2,yr2,xr3,yr3,xr4,yr4,xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,h,k)

R1m_h(k) = sqrt((xr1-xpm_Nh(k))^2+(yr1-ypm_Nh(k))^2+h^2);
R2m_h(k) = sqrt((xr2-xpm_Nh(k))^2+(yr2-ypm_Nh(k))^2+h^2);
R3m_h(k) = sqrt((xr3-xpm_Nh(k))^2+(yr3-ypm_Nh(k))^2+h^2);
R4m_h(k) = sqrt((xr4-xpm_Nh(k))^2+(yr4-ypm_Nh(k))^2+h^2);


r1_partial_x =-1*(xr1-xpm_Nh(k))/R1m_h(k);
r1_partial_y =-1*(yr1-ypm_Nh(k))/R1m_h(k);
r2_partial_x =-1*(xr2-xpm_Nh(k))/R2m_h(k); 
r2_partial_y =-1*(yr2-ypm_Nh(k))/R2m_h(k);
r3_partial_x =-1*(xr3-xpm_Nh(k))/R3m_h(k);
r3_partial_y =-1*(yr3-ypm_Nh(k))/R3m_h(k);
r4_partial_x =-1*(xr4-xpm_Nh(k))/R4m_h(k); 
r4_partial_y =-1*(yr4-ypm_Nh(k))/R4m_h(k);

H = [r1_partial_x 0 r1_partial_y 0 0 0
              r2_partial_x 0 r2_partial_y 0 0 0
              r3_partial_x 0 r3_partial_y 0 0 0
              r4_partial_x 0 r4_partial_y 0 0 0 ];

R = [sig_x_r^2 0 0 0
    0 sig_y_r^2 0 0
    0 0 sig_x_r^2 0
    0 0 0 sig_y_r^2];    
end