function dircos = q2dc(q)

dircos = zeros(3,3);

dircos(1,1) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4);

dircos(2,1) = 2* (q(1)*q(2) - q(3)*q(4));

dircos(3,1) = 2* (q(1)*q(3) + q(2)*q(4));

dircos(1,2) = 2* (q(1)*q(2) + q(3)*q(4));

dircos(2,2) = - q(1)*q(1) + q(2)*q(2) - q(3)*q(3) + q(4)*q(4);

dircos(3,2) = 2* (q(2)*q(3) - q(1)*q(4));

dircos(1,3) = 2* (q(1)*q(3) - q(2)*q(4));

dircos(2,3) = 2* (q(2)*q(3) + q(1)*q(4));

dircos(3,3) = - q(1)*q(1) - q(2)*q(2) + q(3)*q(3) + q(4)*q(4);

