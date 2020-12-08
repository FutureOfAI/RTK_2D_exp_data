function [Vx_h,Vy_h,bx_h,by_h]=upon_optimalflow_measurement(x_update,y_update,k,Vx_h,Vy_h,bx_h,by_h)   
%% 更新資料
Vx_h(k)=Vx_h(k)+x_update(1);
Vy_h(k)=Vy_h(k)+y_update(1);

bx_h(k)=bx_h(k)+x_update(2);
by_h(k)=by_h(k)+y_update(2);

end