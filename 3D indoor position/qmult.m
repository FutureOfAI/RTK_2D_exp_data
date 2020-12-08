function r = qmult(p,q)
% modified from "/home/hs702/thomas/noyola/noyola/qmult.m"
% compute quaternion multification r = p (x) q

qw1 = [ p(4), -p(3),  p(2),  p(1)];
qw2 = [ p(3),  p(4), -p(1),  p(2)];
qw3 = [-p(2),  p(1),  p(4),  p(3)];
qw4 = [-p(1), -p(2), -p(3),  p(4)];
qw  = [qw1;qw2;qw3;qw4];
q   = [ q(1), q(2), q(3), q(4)]';  
r   = qw*q;
