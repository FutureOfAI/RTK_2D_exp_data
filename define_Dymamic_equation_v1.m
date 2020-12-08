function [phi,Q_x,Q_y]=define_Dymamic_equation1(sig_bx,sig_by,sig_xr,sig_yr,dt)
%F = zeros(2);
F = [0 -1
     0  0]; 
%F(2,7)= 0;
%F(5,7)= 0;
phi = expm(F*dt);
q_11x = sig_bx^2;
q_22x = sig_xr^2;
q_11y = sig_by^2;
q_22y = sig_yr^2;

%Q = zeros(2);
Q_x =[(q_11x*dt)+((dt)^3/3)*q_22x -((dt)^2/2)*q_22x
      -((dt)^2/2)*q_22x q_22x*dt];
Q_y =[(q_11y*dt)+((dt)^3/3)*q_22y -((dt)^2/2)*q_22y
      -((dt)^2/2)*q_22y q_22y*dt];
end