function [Rm_data,Rm_Index,nvx_r,nvy_r] = radio_sensor_m_mult(coordinate,x_p_N,y_p_N,sig_x_r,sig_y_r,n,m)
    nvx_r=normrnd(0,sig_x_r,1,n);
    nvy_r=normrnd(0,sig_y_r,1,n);
    Rm=zeros(m);
    [col row] = size(coordinate);
    coor_diff = zeros(col,3);
    for i=1:n
        %% simplify coding
        coor_diff(:,1) = coordinate(:,2)-x_p_N(i);
        coor_diff(:,2) = coordinate(:,3)-y_p_N(i);
        coor_diff(:,3) = coordinate(:,4).*coordinate(:,4);
        Rm(i,1:col) = sqrt(coor_diff(:,1).*coor_diff(:,1) + coor_diff(:,2).*coor_diff(:,2) + coor_diff(:,3)) + nvx_r(i);  
    end
    %% sort most nearest 4 distance
    [Rm_data,Rm_Index] = sort(Rm,2);
end