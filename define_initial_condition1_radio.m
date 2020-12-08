function[x_h,y_h,P00_z]=define_initial_condition1_radio(vxm_h,vym_h,bx_h,by_h,ba0,bp0,g)

%% 宣告狀態變數
%d2r=pi/180;
x_h=zeros(2,1);      % x=[x1;x2;x3;x4;x5;x6;x7;x8]
y_h=zeros(2,1);
%% 定義狀態之初始條件
x_h(1,1) = 0;            % x1(1),velocity error in x-axis (m)
x_h(2,1) = ba0;        % x2(1),x bias
y_h(1,1) = 0;            % x3(1),velocity error in y-axis (m)
y_h(2,1) = bp0;        % x4(1),y bias

x1(1) = x_h(1,1);  
x2(1) = x_h(2,1);  
y1(1) = y_h(1,1);  
y2(1) = y_h(2,1);  

p1 = x1(1)-vxm_h;
p2 = x2(1)-bx_h;
p3 = y1(1)-vym_h;
p4 = y2(1)-by_h;

pii_x = [p1 p2];
pii_x = pii_x.^2;
%P00_x = [pii_x(1) 0
%         0 pii_x(2)];   % diag
P00_x = [pii_x(1) 0
         0 ba0^2];   % diag


pii_y = [p3 p4];
pii_y = pii_y.^2;
%P00_y = [pii_y(1) 0
%         0 pii_y(2)];   % diag
%
P00_y = [pii_y(1) 0
         0 bp0^2];   % diag

end