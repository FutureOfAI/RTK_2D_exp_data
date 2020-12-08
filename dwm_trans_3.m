%
close all;
clear all;
clc;
load('multisapce_data_1070305_v2.mat');%linvelgyro125.mat  linveltest (7)  data (6)aj_1aj_1data (4)
% load('static_test_data_1012.mat');
%aj_7.mat  (y_anchor_2 = 9.37; x_anchor_2 = 0.55;z_anchor_2 = 1.83;)

count = 4000;
count_for = count-1;
data_new = zeros (32,count);
distance_cnt = 1;
inertial_cnt = 1;
t = 1:count;
acc = zeros(3,count);
gro = zeros(3,count);
pos = zeros(4,count);
true_trajectory = 4.64;%6.75;
p =1;% not1 not2 not3 not4 
for j_1=1:count_for
    for i_1=1:32 
        %%
        if data(i_1,j_1) == 11
            if i_1+10 <= 32 && data(i_1+10,j_1) == 128
                if rem(sum( data(i_1:(i_1+8),j_1) ),256) == data(i_1+9,j_1)
                    data_new(1:11,distance_cnt) = data(i_1:(i_1+10),j_1);
                    distance_cnt = distance_cnt+1;
                end
            end
            if i_1+10 < 34&& i_1+10 > 32 && data(i_1-22,j_1+1) == 128
                if  rem(sum( data(i_1:i_1+8,j_1) ) ,256) == data(i_1+9,j_1)                    
                    data_new(1:(33-i_1),distance_cnt) = data(i_1:32,j_1);
                    data_new((34-i_1):11,distance_cnt) = data(1:(i_1-22),j_1+1);
                    distance_cnt = distance_cnt+1;
                end
            end
            if i_1+10 <= 34&& i_1+10 > 33 && data(i_1-22,j_1+1) == 128
                if  rem(sum( data(i_1:i_1+8,j_1) ) ,256) == data(i_1-23,j_1+1)                    
                    data_new(1:(33-i_1),distance_cnt) = data(i_1:32,j_1);
                    data_new((34-i_1):11,distance_cnt) = data(1:(i_1-22),j_1+1);
                    distance_cnt = distance_cnt+1;
                end
            end
            
            if i_1+10 > 34 && data(i_1-22,j_1+1) == 128
                if  rem(sum( data(i_1:32,j_1) ) + sum( data(1:(i_1-24),j_1+1) ),256) == data((i_1-23),j_1+1)
                    data_new(1:(33-i_1),distance_cnt) = data(i_1:32,j_1);
                    data_new((34-i_1):11,distance_cnt) = data(1:(i_1-22),j_1+1);
                    distance_cnt = distance_cnt+1;
                end                
            end
        end
        %%bitand(a, uint16(377));bitand(data(i+19,j),uint16(255))   bitand(data((i-11),j+1),uint16(255))
        if data(i_1,j_1) == 8
            if i_1+20 <= 32 && data(i_1+20,j_1) == 128
                if  rem(sum( data(i_1:(i_1+18),j_1) ),256) == data(i_1+19,j_1);
                    data_new(12:32,inertial_cnt) = data(i_1:(i_1+20),j_1);
                    inertial_cnt = inertial_cnt+1;
                end
            end
            if i_1+20 <= 33 && i_1+20 > 32 && data(i_1-12,j_1+1) == 128
                if  rem(sum( data(i_1:(i_1+18),j_1) ),256) == data(i_1+19,j_1);
                    data_new(12:(44-i_1),inertial_cnt) = data(i_1:32,j_1);
                    data_new((45-i_1):32,inertial_cnt) = data(1:(i_1-12),j_1+1);
                    inertial_cnt = inertial_cnt+1;
                end                
            end       
            
            if i_1+20 <= 34 && i_1+20 > 33 && data(i_1-12,j_1+1) == 128
                if  rem(sum( data(i_1:(i_1+18),j_1) ),256) == data(i_1+13,j_1+1);
                    data_new(12:(44-i_1),inertial_cnt) = data(i_1:32,j_1);
                    data_new((45-i_1):32,inertial_cnt) = data(1:(i_1-12),j_1+1);
                    inertial_cnt = inertial_cnt+1;
                end                
            end       
            if i_1+20 > 34 && data(i_1-12,j_1+1) == 128
                if rem((sum( data(i_1:32,j_1) ) + sum( data(1:(i_1-14),j_1+1))),256)== data((i_1-13),j_1+1)
                    data_new(12:(44-i_1),inertial_cnt) = data(i_1:32,j_1);
                    data_new((45-i_1):32,inertial_cnt) = data(1:(i_1-12),j_1+1);
                    inertial_cnt = inertial_cnt+1;
                end                
            end           
        end
%     if inertial_cnt+1 ~= distance_cnt
%       w=1;  
%     end
    end
end
for j_1=1:count_for
%     for i_1=1:30 
        gro(1,j_1) = (bitshift(data_new(13,j_1),8)+data_new(14,j_1));% gro x
        gro(2,j_1) = (bitshift(data_new(15,j_1),8)+data_new(16,j_1));% gro y
        gro(3,j_1) = (bitshift(data_new(17,j_1),8)+data_new(18,j_1));% gro z
        acc(1,j_1) = (bitshift(data_new(19,j_1),8)+data_new(20,j_1));% acc x
        acc(2,j_1) = (bitshift(data_new(21,j_1),8)+data_new(22,j_1));% acc y
        acc(3,j_1) = (bitshift(data_new(23,j_1),8)+data_new(24,j_1));% acc z  
        pos(1,j_1) = (bitshift(data_new(2,j_1),8)+data_new(3,j_1))/1024;
        pos(2,j_1) = (bitshift(data_new(4,j_1),8)+data_new(5,j_1))/1024;
        pos(3,j_1) = (bitshift(data_new(6,j_1),8)+data_new(7,j_1))/1024;
        pos(4,j_1) = (bitshift(data_new(8,j_1),8)+data_new(9,j_1))/1024;
%     end
end
%% clear bad distance data
for j_1= 2:1:count_for
    
    if(pos(1,j_1)>15)
        pos(1,j_1) = pos(1,j_1-1);
    end
    
    if(pos(2,j_1)>15)
        pos(2,j_1) = pos(2,j_1-1);
    end
    
    if(pos(3,j_1)>15)
        pos(3,j_1) = pos(3,j_1-1);
    end
    
    if(pos(4,j_1)>15)
        pos(4,j_1) = pos(4,j_1-1);
    end
end

for k = 1:count_for        
    for m = 1:3
        if acc(m,k) >= 2^15
            acc(m,k) = acc(m,k) - 2^16;
        end    
        if gro(m,k) >= 2^15
            gro(m,k) = gro(m,k) - 2^16;
        end               
        acc(m,k) = acc(m,k)/8192; %I don't know why  ADC to real?      
        gro(m,k) = (gro(m,k)/16.384*(pi/180));        
%         flow(m,k) = flow(m,k)/bitshift(1,16); %I don't know why  ADC to real?
    end    
end
z_tag =1.2;
%% Least square method
%% aj_9 °ª«×§ïÅÜ
y_anchor_1 = 0.9;
x_anchor_1 = 0;
z_anchor_1 = 2.75;
tri_high_1 = z_anchor_1 - z_tag;
% x_anchor_2 = 9;% x_anchor_2 = 8.91; 
% y_anchor_2 = 0;% y_anchor_2 = 0.06;
% z_anchor_2 = 2.5;% z_anchor_2 = 2.5;


%% aj_9 Àð¾À«e­±
y_anchor_2 = 11.8; 
x_anchor_2 = 0.17;
z_anchor_2 = 2.72;%1.83;
tri_high_2 = z_anchor_2 - z_tag;

y_anchor_3 = 11.7;
x_anchor_3 = 7.85;%7.73 7.78     7.85
z_anchor_3 = 2.5;
tri_high_3 = z_anchor_3 - z_tag;

y_anchor_4 = 0.13;
x_anchor_4 = 7.65;%7.61
z_anchor_4 = 2.7;
tri_high_4 = z_anchor_4 - z_tag;

x_tag = 4.64; y_tag = 6.75; z_tag = 1.32; %x_tag = 4.69
% x_tag = 4.69; y_tag = 6.75; z_tag = 1.32; %x_tag = 4.69
% x_tag = 3.02; y_tag = 5.4; z_tag = 1.2; 
d1 = sqrt((x_tag-x_anchor_1)^2+(y_tag-y_anchor_1)^2+(z_tag-z_anchor_1)^2) ;
d2 = sqrt((x_tag-x_anchor_2)^2+(y_tag-y_anchor_2)^2+(z_tag-z_anchor_2)^2) ;
d3 = sqrt((x_tag-x_anchor_3)^2+(y_tag-y_anchor_3)^2+(z_tag-z_anchor_3)^2) ;
d4 = sqrt((x_tag-x_anchor_4)^2+(y_tag-y_anchor_4)^2+(z_tag-z_anchor_4)^2) ;
% d5 = sqrt((x_tag-x_anchor_5)^2+(y_tag-y_anchor_5)^2) ;

position =zeros(3,(count_for+1)/2);
l=1;
for k = 1:count_for        


distance_sum = -pos(4,k)^2+y_anchor_4^2+x_anchor_4^2+z_anchor_4^2;

b=[pos(1,k)^2-y_anchor_1^2-x_anchor_1^2-z_anchor_1^2+distance_sum
   pos(2,k)^2-y_anchor_2^2-x_anchor_2^2-z_anchor_2^2+distance_sum
   pos(3,k)^2-y_anchor_3^2-x_anchor_3^2-z_anchor_3^2+distance_sum];
%  d4^2-x_anchor_4^2-y_anchor_4^2+distance_sum
A=-2.*[y_anchor_1-y_anchor_4 x_anchor_1-x_anchor_4 z_anchor_1-z_anchor_4 
       y_anchor_2-y_anchor_4 x_anchor_2-x_anchor_4 z_anchor_2-z_anchor_4 
       y_anchor_3-y_anchor_4 x_anchor_3-x_anchor_4 z_anchor_3-z_anchor_4];
%    x_anchor_4-x_anchor_5 y_anchor_4-y_anchor_5


position(:,l)=inv(A'*A)*A'*b;
l=l+1;
end

l = 1;
%% LS
% u=1;
% for  k = 1:count_for    
% H = [(x_anchor_4-x_anchor_1) (y_anchor_4-y_anchor_1) (z_anchor_4-z_anchor_1)
%     (x_anchor_4-x_anchor_2) (y_anchor_4-y_anchor_2) (z_anchor_4-z_anchor_2)
%     (x_anchor_4-x_anchor_3) (y_anchor_4-y_anchor_3) (z_anchor_4-z_anchor_3)];
% r1_n =  pos(1,k)^2 - (x_anchor_1^2+y_anchor_1^2+z_anchor_1^2) ;
% r2_n =  pos(2,k)^2 - (x_anchor_2^2+y_anchor_2^2+z_anchor_2^2) ;
% r3_n =  pos(3,k)^2 - (x_anchor_3^2+y_anchor_3^2+z_anchor_3^2) ;
% r4_n =  pos(4,k)^2 - (x_anchor_4^2+y_anchor_4^2+z_anchor_4^2) ;
% z = 1/2.*[r1_n - r4_n
%          r2_n - r4_n
%          r3_n - r4_n];
% 
% position_4(:,u)=inv(H)*z;
% u=u+1;
% end

%%  LS lterative method
% p_t(1,1,1) = x_tag;
% p_t(2,1,1) = y_tag;
% p_t(3,1,1) = z_tag;
% % r_t(1,1,1) = d1;
% % r_t(2,1,1) = d2;
% % r_t(3,1,1) = d3;
% for k = 1:count_for 
% r_t(1,1,k) = sqrt((p_t(1,1,k)-x_anchor_1)^2+(p_t(2,1,k)-y_anchor_1)^2+(p_t(3,1,k)-z_anchor_1)^2) ;
% r_t(2,1,k) = sqrt((p_t(1,1,k)-x_anchor_2)^2+(p_t(2,1,k)-y_anchor_2)^2+(p_t(3,1,k)-z_anchor_2)^2) ;
% r_t(3,1,k) = sqrt((p_t(1,1,k)-x_anchor_4)^2+(p_t(2,1,k)-y_anchor_4)^2+(p_t(3,1,k)-z_anchor_4)^2) ;
% r_e(1,1,k) = pos(1,k);
% r_e(2,1,k) = pos(2,k);
% r_e(3,1,k) = pos(4,k);
% r_d(1,1,k) = r_e(1,1,k)-r_t(1,1,k);
% r_d(2,1,k) = r_e(2,1,k)-r_t(2,1,k);
% r_d(3,1,k) = r_e(3,1,k)-r_t(3,1,k);
% 
% H_t = [(p_t(1,1,k)-x_anchor_1)/r_e(1,1,k) (p_t(2,1,k)-y_anchor_1)/r_e(1,1,k) (p_t(3,1,k)-z_anchor_1)/r_e(1,1,k)
%     (p_t(1,1,k)-x_anchor_2)/r_e(2,1,k) (p_t(2,1,k)-y_anchor_2)/r_e(2,1,k) (p_t(3,1,k)-z_anchor_2)/r_e(2,1,k)
%     (p_t(1,1,k)-x_anchor_4)/r_e(3,1,k) (p_t(2,1,k)-y_anchor_4)/r_e(3,1,k) (p_t(3,1,k)-z_anchor_4)/r_e(3,1,k)];
% p_t(:,1,k+1) = p_t(:,1,k)+inv(H_t)*r_d(:,1,k);
% p_x(k) =p_t(2,1,k);
% p_y(k) =p_t(1,1,k);
% end

%% Trilateration method
position_2 =zeros(3,(count_for+1)/2);
for k = 1:count_for  
% if(pos(3,k)<pos(1,k) && pos(3,k)<pos(2,k) && pos(3,k)<pos(4,k))
if(p == 3)% not anchor 3
a_1 = y_anchor_4 - y_anchor_1;
b_1 = x_anchor_4 - x_anchor_1;
c_1 = y_anchor_2 - y_anchor_1;
position_2(1,l) = (pos(1,k)^2-pos(2,k)^2+c_1^2) / (2*c_1);
position_2(2,l) = ((pos(1,k)^2-pos(4,k)^2+(position_2(1,l)-a_1)^2)/(2*b_1) + (b_1/2) - ((pos(1,k)^2-pos(2,k)^2+c_1^2)^2)/(8*b_1*(c_1^2)));
position_2(3,l) = sqrt(pos(1,k)^2-position_2(1,l)^2-position_2(2,l)^2);
position_2(3,l) = z_anchor_1 - position_2(3,l);
position_2(1,l) = position_2(1,l) + y_anchor_1;
position_2(2,l) = position_2(2,l) + x_anchor_1;
end
%
% if(pos(4,k)<pos(1,k) && pos(4,k)<pos(2,k) && pos(4,k)<pos(3,k))
if(p == 4)% not anchor 4
% sqrt((y_anchor_3 - y_anchor_1)^2 +(y_anchor_3 - y_anchor_1)^2 + (z_anchor_3 - z_anchor_1)^2)    
%     
a_1 = y_anchor_3 - y_anchor_1;
b_1 = x_anchor_3 - x_anchor_1;
c_1 = y_anchor_2 - y_anchor_1;
position_2(1,l) = (pos(1,k)^2-pos(2,k)^2+c_1^2) / (2*c_1);%+(tri_high_2^2-tri_high_1^2)
position_2(2,l) = ((pos(1,k)^2-pos(3,k)^2+(position_2(1,l)-a_1)^2)/(2*b_1) + (b_1/2) - ((pos(1,k)^2-pos(2,k)^2+c_1^2)^2)/(8*b_1*(c_1^2)));%+(tri_high_2^2-tri_high_1^2)
position_2(3,l) = sqrt(pos(1,k)^2-position_2(1,l)^2-position_2(2,l)^2);
position_2(3,l) = z_anchor_1 - position_2(3,l);
position_2(1,l) = position_2(1,l) + y_anchor_1;
position_2(2,l) = position_2(2,l) + x_anchor_1;
end
%
% if(pos(1,k)<pos(2,k) && pos(1,k)<pos(3,k) && pos(1,k)<pos(4,k))
if(p == 1)% not anchor 1
a_1 = y_anchor_2 - y_anchor_3;
b_1 = x_anchor_2 - x_anchor_3;
c_1 = y_anchor_4 - y_anchor_3;
position_2(1,l) = (pos(3,k)^2-pos(4,k)^2+c_1^2) / (2*c_1);
position_2(2,l) = ((pos(3,k)^2-pos(2,k)^2+(position_2(1,l)-a_1)^2)/(2*b_1) + (b_1/2) - ((pos(3,k)^2-pos(4,k)^2+c_1^2)^2)/(8*b_1*(c_1^2)));
position_2(3,l) = sqrt(pos(3,k)^2-position_2(1,l)^2-position_2(2,l)^2);
position_2(3,l) = z_anchor_3 - position_2(3,l);
position_2(1,l) = position_2(1,l) + y_anchor_3;
position_2(2,l) = position_2(2,l) + x_anchor_3;
end

%
% if(pos(2,k)<pos(1,k) && pos(2,k)<pos(3,k) && pos(2,k)<pos(4,k))
if(p == 2)% not anchor 2
a_1 = y_anchor_1 - y_anchor_3;
b_1 = x_anchor_1 - x_anchor_3;
c_1 = y_anchor_4 - y_anchor_3;
position_2(1,l) = (pos(3,k)^2-pos(4,k)^2+c_1^2) / (2*c_1);
position_2(2,l) = ((pos(3,k)^2-pos(1,k)^2+(position_2(1,l)-a_1)^2)/(2*b_1) + (b_1/2) - ((pos(3,k)^2-pos(4,k)^2+c_1^2)^2)/(8*b_1*(c_1^2)));
position_2(3,l) = sqrt(pos(3,k)^2-position_2(1,l)^2-position_2(2,l)^2);
position_2(3,l) = z_anchor_3 - position_2(3,l);
position_2(1,l) = position_2(1,l) + y_anchor_3;
position_2(2,l) = position_2(2,l) + x_anchor_3;
end
l=l+1;
end

%% Trilateration method 3 Anchor

position_3 =zeros(3,(count_for+1)/2);
l = 1;
for k = 1:count_for  

a_1 = y_anchor_4 - y_anchor_1;
b_1 = x_anchor_4 - x_anchor_1;
c_1 = y_anchor_2 - y_anchor_1;
position_3(1,l) = (pos(1,k)^2-pos(2,k)^2+c_1^2) / (2*c_1);
position_3(2,l) = ((pos(1,k)^2-pos(4,k)^2+(position_3(1,l)-a_1)^2)/(2*b_1) + (b_1/2) - ((pos(1,k)^2-pos(2,k)^2+c_1^2)^2)/(8*b_1*(c_1^2)));
position_3(3,l) = sqrt(pos(1,k)^2-position_3(1,l)^2-position_3(2,l)^2);
position_3(3,l) = z_anchor_1 - position_3(3,l);
position_3(1,l) = position_3(1,l) + y_anchor_1;
position_3(2,l) = position_3(2,l) + x_anchor_1;


l=l+1;
end
%% position differential becone velocity
velocity=zeros(2,count_for);
for i = 2:count_for
    velocity(1,i-1) = (position(1,i) - position(1,i-1))/0.01;% X Positioning to X Velocity / differential
    velocity(2,i-1) = (position(2,i) - position(2,i-1))/0.01;% Y Positioning to Y Velocity / differential
end

%% angular rate integral  become angle
gro_dot=gro';
Euler_angle = zeros(size(gro_dot));
Euler_angle1 = zeros(size(gro_dot));
e = 1;
for i = 2:length(gro_dot)
%     gro_dot(i,:) = (1+0.1909)*(gro_dot(i,:) - (4.37*(pi/180)));
%     gro_dot(i-1,:) = (1+0.1909)*(gro_dot(i-1,:) - (4.37*(pi/180)));
% Euler_angle(i,:) = Euler_angle(i-1,:) + (gro_dot(i,:)) * 0.01;
%     gro_dot_h(i,:) = (1-(-0.0015))*(gro_dot(i,:)-(-4.37*(pi/180)));% °f-0.0512 ¶¶-0.2577
%     gro_dot_h(i-1,:) = (1-(-0.0015))*(gro_dot(i-1,:)-(-4.37*(pi/180)));% °f-0.0512 ¶¶-0.2577
%     Euler_angle(i,:) = Euler_angle(i-1,:)+(gro_dot_h(i,:)+gro_dot_h(i-1,:))*0.01/2;
%% 6/16 compare
%     gro_dot_h(i,:) = (1)*(gro_dot(i,:)-(-4.49*(pi/180)));% °f-0.0512 ¶¶-0.2577
%     gro_dot_h(i-1,:) = (1)*(gro_dot(i-1,:)-(-4.49*(pi/180)));% °f-0.0512 ¶¶-0.2577
%     Euler_angle(i,:) = Euler_angle(i-1,:)+(gro_dot_h(i,:)+gro_dot_h(i-1,:))*0.01/2;
%     
%     gro_dot_h1(i,:) = (1)*(gro_dot(i,:));% °f-0.0512 ¶¶-0.2577
%     gro_dot_h1(i-1,:) = (1)*(gro_dot(i-1,:));% °f-0.0512 ¶¶-0.2577
%     Euler_angle1(i,:) = Euler_angle1(i-1,:)+(gro_dot_h1(i,:)+gro_dot_h1(i-1,:))*0.01/2;
  
%% 6/17 scale factor compare
     gro_dot_h(i,:) = (1-(-0.0015))*(gro_dot(i,:)-(-4.37*(pi/180)));% °f-0.0512 ¶¶-0.2577
     gro_dot_h(i-1,:) = (1-(-0.0015))*(gro_dot(i-1,:)-(-4.37*(pi/180)));% °f-0.0512 ¶¶-0.2577
     Euler_angle(i,:) = Euler_angle(i-1,:)+(gro_dot_h(i,:)+gro_dot_h(i-1,:))*0.01/2;%   

     gro_dot_h1(i,:) = (1)*(gro_dot(i,:)-(-4.37*(pi/180)));% °f-0.0512 ¶¶-0.2577
     gro_dot_h1(i-1,:) = (1)*(gro_dot(i-1,:)-(-4.37*(pi/180)));% °f-0.0512 ¶¶-0.2577
     Euler_angle1(i,:) = Euler_angle1(i-1,:)+(gro_dot_h1(i,:)+gro_dot_h1(i-1,:))*0.01/2;%   
end

%% test acc & gyro (noise,bias)
std_acc_x = std(acc(1,:));%x noise
std_acc_y = std(acc(2,:));%y noise
std_gyro_z = std(gro(3,:));%z noise
% sig_bx_0 = std_x;%err_factor*0.01*g;              % accel noise in g in along-direction
% sig_by_0 = std_y;%err_factor*0.01*g;              % accel noise in g in penpenticular-direction
%2*filtCutOff)/(1/samplePeriod)
mean_acc_x = mean(acc(1,:)); %bias
mean_acc_y = mean(acc(2,:)); %bias
mean_gyro_z = mean(gro(3,:)); %bias

%% low pass filter
samplePeriod = 0.01;  % samplePeriod
order = 1;            % order
filtCutOff =2;       % cutoff frequency
[c, a] = butter(order, filtCutOff/(1/samplePeriod/2), 'low'); %  fc/(fs/2)
low_gro = filtfilt(c, a, gro')';
low_acc = filtfilt(c, a, acc')';
%% plot line
line([1,10],[1.12 1.12])

y = 1.32+zeros(10,1);
y_1 = 1:count;
x_1 = 1.32+zeros(count,1);
x = 1:10;

% x = 1:10;             
% x_1 = 1:count/2;
% y_1 = 1.12+zeros(count/2,1);
% y = 1.2+zeros(10,1);

%% positioning error & mean
for i = count/2:count-1
mean_1(1,i) = true_trajectory-position(2,i);
mean_2(1,i) = true_trajectory-position_2(2,i); 
mean_3(1,i) = 1.32-position_2(3,i); 
end
mean_LS = mean(mean_1(1,count/2:count-1));
mean_Tri = mean(mean_2(1,count/2:count-1));
mean_Tri_h = mean(mean_3(1,count/2:count-1));
std_LS = std(mean_1(1,count/2:count-1));
std_Tri = std(mean_2(1,count/2:count-1));
std_Tri_1 = std(mean_3(1,count/2:count-1));

% bias_err+3*std
sigma_LS = abs(mean_LS) + (3*std_LS);
sigma_Tri = abs(mean_Tri) + (3*std_Tri);
max_1(1,1) = max(abs(mean_1(1,:)));
max_1(1,2) = min(abs(mean_1(1,:)));
%% distance error & mean
for i = 1:count
mean_d_(1,i) = d1-pos(1,i);
mean_d_(2,i) = d2-pos(2,i); 
mean_d_(3,i) = d3-pos(3,i);
mean_d_(4,i) = d4-pos(4,i); 
end
mean_d(1,1) = mean(mean_d_(1,:));
mean_d(2,1) = mean(mean_d_(2,:));
mean_d(3,1) = mean(mean_d_(3,:));
mean_d(4,1) = mean(mean_d_(4,:));

std_d(1,1) = std(mean_d_(1,:));
std_d(2,1) = std(mean_d_(2,:));
std_d(3,1) = std(mean_d_(3,:));
std_d(4,1) = std(mean_d_(4,:));
% bias_err+3*std
sigma_d(1,1) = abs(mean_d(1,1)) + (3*std_d(1,1));
sigma_d(2,1) = abs(mean_d(2,1)) + (3*std_d(2,1));
sigma_d(3,1) = abs(mean_d(3,1)) + (3*std_d(3,1));
sigma_d(4,1) = abs(mean_d(4,1)) + (3*std_d(4,1));

%% CDF (Cumulative distribution functions)

pd_Tri = makedist('Normal',mean_Tri,std_Tri);%('Normal',mean,std)
pd_LS = makedist('Normal',mean_LS,std_LS);%('Normal',mean,std)
ranger = 0:0.01:0.35;
cdf_Tri = cdf(pd_Tri,ranger);
cdf_LS = cdf(pd_LS,ranger);


%% circle

r_circle = 1.24; % ¶ê¥b®|
x_circle = 4.17; % ¶ê¤ß x ®y¼Ð
y_circle = 2.46; % ¶ê¤ß y ®y¼Ð
theta=0:pi/50:2*pi;
x_circle_1=x_circle+r_circle*cos(theta);
y_circle_1=y_circle+r_circle*sin(theta);
%%
figure (1)
plot(y,x,'--r',position(1,:),position(2,:),y_anchor_1,x_anchor_1,'r*')
hold on
% plot(5.4,1.4,'diamond',5.4,9.2,'rsquare',5.4,1.4,'bx','MarkerSize',13,'LineWidth',2) % 10 ¨Ó¦^
% plot(5.4,0.74,'diamond',5.4,6.8,'rsquare','MarkerSize',13,'LineWidth',2)%,5.4,0.74,'bx'
plot(0.85,3.35,'diamond',5.45,4.11,'rsquare','MarkerSize',13,'LineWidth',2)%,5.4,0.74,'bx'
legend('Reference path','Least square','Anchor','Start point','End point');%'Reentry point',
hold on
plot(y_anchor_2,x_anchor_2,'r*',y_anchor_3,x_anchor_3,'r*',y_anchor_4,x_anchor_4,'r*')
title('Least square')
xlabel('Y position in m')
ylabel('X position in m')
axis([-2 13 -2 10])
set(gca,'YDir','reverse')
grid

figure (2)
plot(y,x,'--r',position_2(1,:),position_2(2,:),y_anchor_1,x_anchor_1,'r*')
hold on
plot(5.4,0.74,'diamond',5.4,6.8,'rsquare',5.4,0.74,'bx','MarkerSize',13,'LineWidth',2)
legend('Reference path','Least square','Anchor','Start point','Reentry point','End point');
hold on
plot(y_anchor_2,x_anchor_2,'r*',y_anchor_3,x_anchor_3,'r*',y_anchor_4,x_anchor_4,'r*')
anchor_d=p;
title(['Trilateration method not Anchor',num2str(anchor_d)])
xlabel('Y position in m')
ylabel('X position in m')
axis([-2 13 -2 10])
set(gca,'YDir','reverse')
grid
figure (3)
plot(y,x,'--r',position_3(1,:),position_3(2,:),y_anchor_1,x_anchor_1,'r*')
hold on
plot(5.4,0.74,'diamond',5.4,6.8,'rsquare',5.4,0.74,'bx','MarkerSize',13,'LineWidth',2)
legend('Reference path','Least square','Anchor','Start point','Reentry point','End point');
hold on
plot(y_anchor_2,x_anchor_2,'r*',y_anchor_3,x_anchor_3,'r*',y_anchor_4,x_anchor_4,'r*')
title('Trilateration method 3 Anchor')
xlabel('Y position in m')
ylabel('X position in m')
axis([-2 13 -2 10])
set(gca,'YDir','reverse')
grid
g=9.8;
figure (4)
subplot(211)
plot(t/100,acc(1,:)/g,'b',t/100,low_acc(1,:)/g,'r')
xlabel('Time(s)')
ylabel('Y accelerometer in g')
subplot(212)
plot(t/100,acc(2,:)/g,'r',t/100,low_acc(2,:)/g,'b')
xlabel('Time(s)')
ylabel('X accelerometer in g')


figure (5)
subplot(211)
plot(t/100,gro(3,:)*180/pi,'b')%t,gro(1,:),'b',t,gro(2,:),'r',t,    ,t/100,linVelLP(:,3)*180/pi,'r'
xlabel('Time(s)')
ylabel('Z Gyro in degree/s')
subplot(212)
plot(t/100,(gro(3,:)+0.0824)*180/pi,'b')%t,gro(1,:),'b',t,gro(2,:),'r',t,    ,t/100,linVelLP(:,3)*180/pi,'r'
xlabel('Time(s)')
ylabel('Z Gyro in degree/s')


figure (7)
title('Trilateration method')
subplot(211)
plot(1:count-1,position_2(1,:))
xlabel('time')
ylabel('Y Positioning')
subplot(212)
% plot(x_1,y_1,'r',1:count/2,position_2(2,:))%,1:1500,y__1
plot(y_1,x_1,'r',1:count-1,position_2(2,:))%,1:1500,y__1
xlabel('time')
ylabel('X Positioning')

figure (8)
title('Least square method')
subplot(211)
% plot(1:count/2,position(1,:))
plot(1:count-1,position(1,:))%,1:1500,y__1
xlabel('time')
ylabel('Y Positioning')
subplot(212)
% plot(x_1,y_1,'r',1:count/2,position(2,:))%,1:1500,y__1
plot(y_1,x_1,'r',1:count-1,position(2,:))
xlabel('time')
ylabel('X Positioning')

figure (9)
title('Least square method')
subplot(211)
plot(1:count_for,velocity(1,:))
xlabel('time')
ylabel('Y velocity m/s')
subplot(212)
plot(1:count_for,velocity(2,:))%,1:1500,y__1
xlabel('time')
ylabel('X velocity m/s')
% xlswrite('data1.xls',gro(1,:))

figure (10)
plot(t/100,(Euler_angle(:,3)*180/pi),'b',t/100,(Euler_angle1(:,3)*180/pi),'r');
xlabel('Time(s)');
ylabel('Angle');
grid
%% 6/17  compare
% subplot(211)
% plot(t/100,(Euler_angle1(:,3)*180/pi), 'b');
% xlabel('Time(s)');
% ylabel('Angle');
% grid
% subplot(212)
% plot(t/100,(Euler_angle(:,3)*180/pi), 'b');
% xlabel('Time(s)');
% ylabel('Angle');
% grid



figure (11)
title('Range measurements')
subplot(411)
plot(1:count_for,pos(1,1:count_for))
xlabel('time')
ylabel('Anchor1 distance (m)')
subplot(412)
plot(1:count_for,pos(2,1:count_for))
xlabel('time')
ylabel('Anchor2 distance (m)')
subplot(413)
plot(1:count_for,pos(3,1:count_for))
xlabel('time')
ylabel('Anchor3 distance (m)')
subplot(414)
plot(1:count_for,pos(4,1:count_for))
xlabel('time')
ylabel('Anchor4 distance (m)')

figure (12)
title('axm aym wzm measurements')
subplot(311)
plot(t/100,acc(1,:),'b')
axis([0 count/100 -0.1 0.1])
xlabel('Time(s)')
ylabel('aym in g')
subplot(312)
plot(t/100,acc(2,:),'b')
axis([0 count/100 -0.1 0.1])
xlabel('Time(s)')
ylabel('axm in g')
subplot(313)
plot(t/100,gro(3,:)*180/pi,'b')%t,gro(1,:),'b',t,gro(2,:),'r',t,    ,t/100,linVelLP(:,3)*180/pi,'r'
axis([0 count/100 -10 0])
xlabel('Time(s)')
ylabel('wzm in degree/s')

figure(13)
title('Least square method error')
plot(1:count-1,true_trajectory-position(1,:))%,1:1500,y__1
xlabel('time')
ylabel('Y Positioning error')
%% CDF
figure(14)
title('Cumulative distribution functions')
plot(ranger,cdf_Tri,ranger,cdf_LS)%,1:1500,y__1
legend('Tri','LS')
xlabel('position error[m]')
ylabel('CDF')
grid

%% Circle
figure(15)
title('Circle trajectory')
plot(x_circle_1,y_circle_1)%,1:1500,y__1

xlabel('x[m]')
ylabel('y[m]')
grid



figure (16)
title('Least square method')

plot(1:count-1,position(3,:),'b',y_1,x_1,'r')
xlabel('time')
ylabel('Z Positioning')

%% LS lterative method
% figure (17)
% plot(p_x(1:3999),p_y(1:3999))
% hold on
% title('Least square')
% xlabel('Y position in m')
% ylabel('X position in m')
% axis([-2 13 -2 10])
% set(gca,'YDir','reverse')
% grid

% d1
% d2
% d3
% d4
mean_Tri
std_Tri
sigma_Tri