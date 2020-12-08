function [x_h,y_h,P00_x,P00_y,K_x,K_y,x_update,y_update]=Kalman_Filter_update1(x_h,y_h,P00_x,P00_y,Vx_h,Vy_h,vxm,vym,H_x,H_y,R_x,R_y,k)
%% ( k| k-1) ->( k| k)¡G­×¥¿ª¬ºA¥H¤ÎP
% x(:,k)=[x1(k);...;x8(k)],size of x = 8*n
zxm_x(k) = vxm(k)-Vx_h(k) ;
Mu_x=zxm_x(k);
K_x=P00_x*(H_x')*(inv(H_x*P00_x*(H_x')+R_x));      % P(k|k-1)
%K(7) = 0;
%K(8) = 0;
x_update = K_x*Mu_x;
x_h=x_h+x_update;                        % x(k|k)
I=eye(2);
P00_x=(I-K_x*H_x)*P00_x;                     % P(k|k)


zxm_y(k) = vym(k)-Vy_h(k) ;
Mu_y=zxm_y(k);
K_y=P00_y*(H_y')*(inv(H_y*P00_y*(H_y')+R_y));      % P(k|k-1)
%K(7) = 0;
%K(8) = 0;
y_update = K_y*Mu_y;
y_h=y_h+y_update;                        % x(k|k)
I=eye(2);
P00_y=(I-K_y*H_y)*P00_y;                     % P(k|k)
end

