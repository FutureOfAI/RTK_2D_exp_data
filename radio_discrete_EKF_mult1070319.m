function[H,R,Rm_h]=radio_discrete_EKF_mult1070319(coordinate,Rm_Index,xpm_Nh,ypm_Nh,sig_x_r,sig_y_r,k,format_1)
if (format_1 == 1) % straight line
    
    coor_diff(:,1) = coordinate(Rm_Index(k,1:4),2)-xpm_Nh(k);
    coor_diff(:,2) = coordinate(Rm_Index(k,1:4),3)-ypm_Nh(k);
    coor_diff(:,3) = coordinate(Rm_Index(k,1:4),4)-1.32;
    Rm_h(k,1:4) = sqrt( coor_diff(:,1).*coor_diff(:,1) + coor_diff(:,2).*coor_diff(:,2) + coor_diff(:,3).*coor_diff(:,3) );
    
    r_partial_x = -( coordinate(Rm_Index(k,1:4),2)-xpm_Nh(k) )./Rm_h(k,1:4)';
    r_partial_y = -( coordinate(Rm_Index(k,1:4),3)-ypm_Nh(k) )./Rm_h(k,1:4)';

    H = [r_partial_x zeros(4,1) zeros(4,1) r_partial_y zeros(4,1) zeros(4,1) zeros(4,1) zeros(4,1)];    
elseif (format_1 == 2)  
    R1m_h(k) = sqrt((coordinate(Rm_Index(k,1),2)-xpm_Nh(k))^2+(coordinate(Rm_Index(k,1),3)-ypm_Nh(k))^2+(coordinate(Rm_Index(k,1),4)-1.32)^2); 
    R2m_h(k) = sqrt((coordinate(Rm_Index(k,2),2)-xpm_Nh(k))^2+(coordinate(Rm_Index(k,2),3)-ypm_Nh(k))^2+(coordinate(Rm_Index(k,2),4)-1.32)^2);
    R3m_h(k) = sqrt((coordinate(Rm_Index(k,3),2)-xpm_Nh(k))^2+(coordinate(Rm_Index(k,3),3)-ypm_Nh(k))^2+(coordinate(Rm_Index(k,3),4)-1.32)^2);
    R4m_h(k) = sqrt((coordinate(Rm_Index(k,4),2)-xpm_Nh(k))^2+(coordinate(Rm_Index(k,4),3)-ypm_Nh(k))^2+(coordinate(Rm_Index(k,4),4)-1.32)^2);
    
    Rm_h(k,1:4) = [R1m_h(k);R2m_h(k);R3m_h(k);R4m_h(k)];
    
    r1_partial_x =-(coordinate(Rm_Index(k,1),2)-xpm_Nh(k))/R1m_h(k);
    r1_partial_y =-(coordinate(Rm_Index(k,1),3)-ypm_Nh(k))/R1m_h(k);
    r2_partial_x =-(coordinate(Rm_Index(k,2),2)-xpm_Nh(k))/R2m_h(k); 
    r2_partial_y =-(coordinate(Rm_Index(k,2),3)-ypm_Nh(k))/R2m_h(k);
    r3_partial_x =-(coordinate(Rm_Index(k,3),2)-xpm_Nh(k))/R3m_h(k);
    r3_partial_y =-(coordinate(Rm_Index(k,3),3)-ypm_Nh(k))/R3m_h(k);
    r4_partial_x =-(coordinate(Rm_Index(k,4),2)-xpm_Nh(k))/R4m_h(k); 
    r4_partial_y =-(coordinate(Rm_Index(k,4),3)-ypm_Nh(k))/R4m_h(k);

    H = [r1_partial_x 0 0 r1_partial_y 0 0 0 0
                  r2_partial_x 0 0 r2_partial_y 0 0 0 0
                  r3_partial_x 0 0 r3_partial_y 0 0 0 0
                  r4_partial_x 0 0 r4_partial_y 0 0 0 0 ];    
end

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