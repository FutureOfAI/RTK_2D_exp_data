function [Tx,Ty] = triangle(R1m,R2m,xr1,yr1,xr2,yr2,h,Yn_1,Xn_1)

X = (R1m^2-R2m^2+h^2)/2*h;% c is anchor & anchor x axis distance 
Y_1 = sqrt(R1m^2-X^2);%positive
Y_2 = -sqrt(R1m^2-X^2);%nagative

a = (Yn_1-yr1)/(Xn_1-xr1) - (yr2-yr1)/(xr2-xr1);%previous X,Y
a1 = (Y_1-yr1)/(X-xr1) - (yr2-yr1)/(xr2-xr1);
% a2 = (X-dsr1_y)/(Y_2-dsr1_y)-(dsr2_y-dsr1_y)/(dsr2_x-dsr1_x);
    if(a*a1>0)
        Tx = xr1 + X;%
        Ty = xr1 + Y_1;
    else
        Tx = xr1 + X;%
        Ty = xr1 + Y_2;        
    end
    
end


