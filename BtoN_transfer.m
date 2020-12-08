function[x_v_N,y_v_N,x_p_N,y_p_N] = BtoN_transfer(x_v,y_v,x_p,y_p,psi,n)


x_p_N = zeros(1,n);
y_p_N = zeros(1,n);

x_v_N = zeros(1,n);
y_v_N = zeros(1,n);

for i = 1:n
    x_v_N(i) = [cos(psi(i)) -sin(psi(i))]*[x_v(i);y_v(i)];
    y_v_N(i) = [sin(psi(i)) cos(psi(i))]*[x_v(i);y_v(i)];


%    x_p_N(1,i) = x_p_N(1,i-1) + x_v_N(i) * dt;% velocity to position / integral (prediction)
%    y_p_N(1,i) = y_p_N(1,i-1) + y_v_N(i) * dt;
    
    x_p_N(i) = [cos(psi(i)) -sin(psi(i))]*[x_p(i);y_p(i)];
    y_p_N(i) = [sin(psi(i)) cos(psi(i))]*[x_p(i);y_p(i)];

    
end

end