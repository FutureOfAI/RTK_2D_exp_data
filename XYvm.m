function [vxm,vym,nvx,nvy]=XYvm(x_v,y_v,m,n,sig_x,sig_y)
nvx=normrnd(0,sig_x,1,n);
nvy=normrnd(0,sig_y,1,n);
vxm=zeros(m);
vym=zeros(m);
for i=1:n
    vxm(i)=x_v(i)+nvx(i);
    vym(i)=y_v(i)+nvy(i);
end
end