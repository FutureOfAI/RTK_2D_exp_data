function [Rm_data_valid,Rm_data_valid_inx] = radio_sensor_m_mult1070415(pos,cell_anchors,n)
    % sort distance data by ascend
%     [Rm_data,Rm_Index] = sort(pos,1);
    [col,row] = size(pos);
    cell = zeros(3,n); % select anchors cell
    dis_tag_anchor = zeros(3,n); % distance between tag and anchors
    % inteligent select anchors
    for i=1:n
        for j=1:col
            if 0<j && j<=4  % define 1~4 cell_1
                if pos(j,i)>0
                    cell(1,i) = cell(1,i)+1;
                    dis_tag_anchor(1,i) = dis_tag_anchor(1,i) + pos(j,i);
                end
            end     
            if 2<j && j<=6 % define 3~6 cell_2
                if pos(j,i)>0
                    cell(2,i) = cell(2,i)+1;
                    dis_tag_anchor(2,i) = dis_tag_anchor(2,i) + pos(j,i);
                end
            end      
            if 6<j && j<=10 % define 7~10 cell_3
                if pos(j,i)>0
                    cell(3,i) = cell(3,i)+1;
                    dis_tag_anchor(3,i) = dis_tag_anchor(3,i) + pos(j,i);
                end
            end                
        end
        
        [dis_data,dis_Index] = sort(dis_tag_anchor(:,i),1);
        dis_data_save(:,i) = dis_data;
        dis_Index_save(:,i) = dis_Index;
        if cell(dis_Index(1),i) == 4
            Rm_data_valid(:,i) = pos(cell_anchors(dis_Index(1),:),i);
            Rm_data_valid_inx(:,i) = cell_anchors(dis_Index(1),:);             
        else if cell(dis_Index(2),i) == 4
                Rm_data_valid(:,i) = pos(cell_anchors(dis_Index(2),:),i);
                Rm_data_valid_inx(:,i) = cell_anchors(dis_Index(2),:);               
            else if cell(dis_Index(3),i) == 4
                    Rm_data_valid(:,i) = pos(cell_anchors(dis_Index(3),:),i);
                    Rm_data_valid_inx(:,i) = cell_anchors(dis_Index(3),:);                       
                else
                    Rm_data_valid(:,i) = Rm_data_valid(:,i-1);
                    Rm_data_valid_inx(:,i) = Rm_data_valid_inx(:,i-1);                        
%                     Rm_data_valid(:,i) = [0, 0 ,0 , 0];
%                     Rm_data_valid_inx(:,i) = [0, 0 ,0 , 0];                        

                end
            end
        end
        
        
%         if length( find(cell(:,i))==4 )==1
%             if cell(1,i) == 4
%                 Rm_data_valid(:,i) = [pos(cell_anchors(1),i),pos(cell_anchors(2),i),pos(cell_anchors(3),i),pos(cell_anchors(4),i)];
%                 Rm_data_valid_inx(:,i) = [cell_anchors(1),cell_anchors(2),cell_anchors(3),cell_anchors(4)];   
%             else if cell(2,i) == 4
%                 Rm_data_valid(:,i) = [pos(cell_anchors(3),i),pos(cell_anchors(4),i),pos(cell_anchors(5),i),pos(cell_anchors(6),i)];
%                 Rm_data_valid_inx(:,i) = [cell_anchors(3),cell_anchors(4),cell_anchors(5),cell_anchors(6)];   
%                 else
%                 Rm_data_valid(:,i) = [pos(cell_anchors(7),i),pos(cell_anchors(8),i),pos(cell_anchors(9),i),pos(cell_anchors(10),i)];
%                 Rm_data_valid_inx(:,i) = [cell_anchors(7),cell_anchors(8),cell_anchors(9),cell_anchors(10)];  
%                 end
%             end
%         else if length( find(cell(:,i))==4 ) > 1
%                 
%             else
%                 Rm_data_valid(:,i) = Rm_data_valid(:,i-1);
%                 Rm_data_valid_inx(:,i) = Rm_data_valid_inx(:,i-1);    
%             end
%         end
    end
    % find valid data
%     for i = 1:n
%        Rm_data_non_zero = find(Rm_data(:,i));
%        if size(Rm_data_non_zero,1)>3
%            Rm_data_valid(:,i) = [Rm_data(Rm_data_non_zero(1),i),Rm_data(Rm_data_non_zero(2),i),Rm_data(Rm_data_non_zero(3),i),Rm_data(Rm_data_non_zero(4),i)];
%            Rm_data_valid_inx(:,i) = [Rm_Index(Rm_data_non_zero(1),i),Rm_Index(Rm_data_non_zero(2),i),Rm_Index(Rm_data_non_zero(3),i),Rm_Index(Rm_data_non_zero(4),i)];
%        else
%            Rm_data_valid(:,i) = Rm_data_valid(:,i-1);
%            Rm_data_valid_inx(:,i) = Rm_data_valid_inx(:,i-1);
%        end
%     end
end