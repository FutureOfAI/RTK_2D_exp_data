function [x_h,y_h,P00_x,P00_y]=define_initial_condition1(vxerr,vyerr,ba0,bp0)
%% 宣告狀態變數
%d2r=pi/180;
x_h=zeros(2,1);      % x=[x1;x2;x3;x4;x5;x6;x7;x8]
y_h=zeros(2,1);
%% 定義狀態之初始條件
x_h(1,1) = vxerr;            % x1(1),velocity error in x-axis (m)
x_h(2,1) = ba0;        % x2(1),x bias
y_h(1,1) = vyerr;            % x3(1),velocity error in y-axis (m)
y_h(2,1) = bp0;        % x4(1),y bias

P00_x = [vxerr^2 0
         0 ba0^2];   % diag

P00_y = [vyerr^2 0
         0 bp0^2];   % diag

end