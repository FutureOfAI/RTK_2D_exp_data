function q = DTQ(b)
% Input: b (3,3), coordinate transformation matrix
% Output: q (4,1) Euler symmetric parameters (quaternions)
% Ref: Wertz, page 414-415
q=zeros(4,1);
q(4) = .5 * sqrt( 1+b(1,1)+b(2,2)+b(3,3) );
q(1) = .5 * sqrt( 1+b(1,1)-b(2,2)-b(3,3) );
q(2) = .5 * sqrt( 1-b(1,1)+b(2,2)-b(3,3) );
q(3) = .5 * sqrt( 1-b(1,1)-b(2,2)+b(3,3) );
[qtemp,i] = sort(q);
if i(4) == 1,
       q(2) = .25 * ( b(2,1) + b(1,2) ) / q(1);
       q(3) = .25 * ( b(3,1) + b(1,3) ) / q(1);
       q(4) = .25 * ( b(2,3) - b(3,2) ) / q(1);
end;
if i(4) == 2,
       q(1) = .25 * ( b(2,1) + b(1,2) ) / q(2);
       q(3) = .25 * ( b(2,3) + b(3,2) ) / q(2);
       q(4) = .25 * ( b(3,1) - b(1,3) ) / q(2);
end;
if i(4) == 3,
       q(1) = .25 * ( b(1,3) + b(3,1) ) / q(3);
       q(2) = .25 * ( b(2,3) + b(3,2) ) / q(3);
       q(4) = .25 * ( b(1,2) - b(2,1) ) / q(3);
end;
if i(4) == 4,
       q(1) = .25 * ( b(2,3) - b(3,2) ) / q(4);
       q(2) = .25 * ( b(3,1) - b(1,3) ) / q(4);
       q(3) = .25 * ( b(1,2) - b(2,1) ) / q(4);
end;
if q(4) < 0,
   q = -q ;
end;
return;
