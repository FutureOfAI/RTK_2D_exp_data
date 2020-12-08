function [thetam]=Thetam1(theta,wzm,m,n,dt)
thetam=zeros(m);
thetam(1)=theta(1);
for i=1:(n-1)
    thetam(i+1)=thetam(i)+wzm(i)*dt;
end
end