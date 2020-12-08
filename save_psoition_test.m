close all;
clear all;
clc;
load('position1.mat');%linvelgyro125.mat  linveltest (7)  data (6)aj_1aj_1data (4)
%aj_7.mat  (y_anchor_2 = 9.37; x_anchor_2 = 0.55;z_anchor_2 = 1.83;)

count = 200;
count_for = count-1;
data_new = zeros (9,count);

t = 1:count;

pos = zeros(3,count);

for j_1=1:count_for
%     for i_1=1:30 
        pos(1,j_1) = (bitshift(data(2,j_1),8)+data(3,j_1))/1024;
        pos(2,j_1) = (bitshift(data(4,j_1),8)+data(5,j_1))/1024;
        pos(3,j_1) = (bitshift(data(6,j_1),8)+data(7,j_1))/1024;

%     end
end


