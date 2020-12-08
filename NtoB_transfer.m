function[x_v_B,y_v_B,x_p_B,y_p_B,x_a_B,y_a_B] = NtoB_transfer(x_p_N,x_v_N,x_a_N,y_p_N,y_v_N,y_a_N,psi,n)


x_p_B = zeros(1,n);
y_p_B = zeros(1,n);

x_v_B = zeros(1,n);
y_v_B = zeros(1,n);

x_a_B = zeros(1,n);
y_a_B = zeros(1,n);

for i = 1:n
    x_v_B(i) = [cos(psi(i)) sin(psi(i))]*[x_v_N(i);y_v_N(i)];
    y_v_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_v_N(i);y_v_N(i)];

    
    x_p_B(i) = [cos(psi(i)) sin(psi(i))]*[x_p_N(i);y_p_N(i)];
    y_p_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_p_N(i);y_p_N(i)];

%    x_a_B(i) = [cos(psi(i)) sin(psi(i))]*[x_a_N(i);y_a_N(i)];
%    y_a_B(i) = [-sin(psi(i)) cos(psi(i))]*[x_a_N(i);y_a_N(i)];
    x_a_B(i) = cos(psi(i))*x_a_N(i) + sin(psi(i))*y_a_N(i);
    y_a_B(i) = -sin(psi(i))*x_a_N(i) + cos(psi(i))*y_a_N(i);
end

end