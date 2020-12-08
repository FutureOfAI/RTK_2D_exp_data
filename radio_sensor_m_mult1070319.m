function [Rm_data_valid,Rm_data_valid_inx] = radio_sensor_m_mult1070319(pos,n)
    % sort distance data by ascend
    [Rm_data,Rm_Index] = sort(pos,1);   
    % find valid data
    for i = 1:n
       Rm_data_non_zero = find(Rm_data(:,i));
       if size(Rm_data_non_zero,1)>3
           Rm_data_valid(:,i) = [Rm_data(Rm_data_non_zero(1),i),Rm_data(Rm_data_non_zero(2),i),Rm_data(Rm_data_non_zero(3),i),Rm_data(Rm_data_non_zero(4),i)];
           Rm_data_valid_inx(:,i) = [Rm_Index(Rm_data_non_zero(1),i),Rm_Index(Rm_data_non_zero(2),i),Rm_Index(Rm_data_non_zero(3),i),Rm_Index(Rm_data_non_zero(4),i)];
       else
           Rm_data_valid(:,i) = Rm_data_valid(:,i-1);
           Rm_data_valid_inx(:,i) = Rm_data_valid_inx(:,i-1);
       end
    end
end