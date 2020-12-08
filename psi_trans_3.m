%%
close all;
clear all;
clc;
load('aj_psi.mat');%linvelgyro125.mat  linveltest (7)  data (6)aj_1aj_1data (4)
%aj_7.mat  (y_anchor_2 = 9.37; x_anchor_2 = 0.55;z_anchor_2 = 1.83;)

count = 2000;
count_for = count-1;
t = 1:count-1;
% data_new = zeros (32,count);
% for j_1=1:count_for
%     for i_1=1:5 
%         if data(i_1,j_1) == 17
%             if i_1+10 <= 32 && data(i_1+10,j_1) == 128
%                 if rem(sum( data(i_1:(i_1+8),j_1) ),256) == data(i_1+9,j_1)
%                     data_new(1:11,distance_cnt) = data(i_1:(i_1+10),j_1);
%                     distance_cnt = distance_cnt+1;
%                 end
%             end
%             if i_1+10 < 34&& i_1+10 > 32 && data(i_1-22,j_1+1) == 128
%                 if  rem(sum( data(i_1:i_1+8,j_1) ) ,256) == data(i_1+9,j_1)                    
%                     data_new(1:(33-i_1),distance_cnt) = data(i_1:32,j_1);
%                     data_new((34-i_1):11,distance_cnt) = data(1:(i_1-22),j_1+1);
%                     distance_cnt = distance_cnt+1;
%                 end
%             end
%             if i_1+10 <= 34&& i_1+10 > 33 && data(i_1-22,j_1+1) == 128
%                 if  rem(sum( data(i_1:i_1+8,j_1) ) ,256) == data(i_1-23,j_1+1)                    
%                     data_new(1:(33-i_1),distance_cnt) = data(i_1:32,j_1);
%                     data_new((34-i_1):11,distance_cnt) = data(1:(i_1-22),j_1+1);
%                     distance_cnt = distance_cnt+1;
%                 end
%             end
% 
%             if i_1+10 > 34 && data(i_1-22,j_1+1) == 128
%                 if  rem(sum( data(i_1:32,j_1) ) + sum( data(1:(i_1-24),j_1+1) ),256) == data((i_1-23),j_1+1)
%                     data_new(1:(33-i_1),distance_cnt) = data(i_1:32,j_1);
%                     data_new((34-i_1):11,distance_cnt) = data(1:(i_1-22),j_1+1);
%                     distance_cnt = distance_cnt+1;
%                 end                
%             end
%         end
%     end
% end

for j_1=1:count_for

        psi(1,j_1) =( bitshift(data(4,j_1),8)+data(5,j_1));% 

end

for k = 1:count_for        

        if psi(k) >= 2^15
            psi(k) = psi(k) - 2^16;
        end    
        psi(k) = (psi(k)/1024); %I don't know why  ADC to real?   *180/pi   
%         flow(m,k) = flow(m,k)/bitshift(1,16); %I don't know why  ADC to real?
 
end


figure (10)


plot(t/100,(psi(:)),'b');
xlabel('sample');
ylabel('Angle');
grid
%text('t','right')

title('Euler angle from gyroscope');
