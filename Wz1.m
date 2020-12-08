function [wz]=Wz1(x2,x4,m,n)
wz=zeros(m);
for i=1:n
    wz(i)=x2(i)+x4(i);   
end
end