function[H,R,Rm_h]=radio_discrete_EKF_mult(coordinate,Rm_Index,xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,k,format_1)
if (format_1 == 1) % straight line
    
    coor_diff(:,1) = coordinate(Rm_Index(k,1:4),2)-xpm_Nh(k);
    coor_diff(:,2) = coordinate(Rm_Index(k,1:4),3)-ypm_Nh(k);
    coor_diff(:,3) = coordinate(Rm_Index(k,1:4),4);
    Rm_h(k,1:4) = sqrt( coor_diff(:,1).*coor_diff(:,1) + coor_diff(:,2).*coor_diff(:,2) + coor_diff(:,3).*coor_diff(:,3) );

elseif (format_1 == 2)% square
%     R1m_h(k) = sqrt((xr1-xpm_Nh(k))^2+(yr1-ypm_Nh(k))^2+(h_1-1.2^2)); 
%     R2m_h(k) = sqrt((xr2-xpm_Nh(k))^2+(yr2-ypm_Nh(k))^2+(h_2-1.2)^2);
%     R3m_h(k) = sqrt((xr3-xpm_Nh(k))^2+(yr3-ypm_Nh(k))^2+(h_3-1.2)^2);
%     R4m_h(k) = sqrt((xr4-xpm_Nh(k))^2+(yr4-ypm_Nh(k))^2+(h_4-1.2)^2);
end

r_partial_x = -( coordinate(Rm_Index(k,1:4),2)-xpm_Nh(k) )./Rm_h(k,1:4)';
r_partial_y = -( coordinate(Rm_Index(k,1:4),3)-ypm_Nh(k) )./Rm_h(k,1:4)';

H = [r_partial_x zeros(4,1) zeros(4,1) r_partial_y zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1)];

% if k < 1500
R = 1*[sig_x_r^2 0 0 0
    0 sig_y_r^2 0 0
    0 0 sig_x_r^2 0
    0 0 0 sig_y_r^2];   
% else
% R = 100*[sig_x_r^2 0 0 0
%         0 sig_y_r^2 0 0
%         0 0 sig_x_r^2 0
%         0 0 0 sig_y_r^2];  
% end

end