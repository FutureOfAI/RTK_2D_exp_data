bias_c = 4.37;
aj_16_n1 = 1-(69.48 + (bias_c*5))/90;
aj_16_n2 = 1-(67.85 + (bias_c*5))/90;
aj_16_n3 = 1-(67.82 + (bias_c*5))/90;
aj_16_n4 = 1-(67.19 + (bias_c*5))/90;
aj_16_n5 = 1-(68.3 + (bias_c*5))/90;


aj_16_p1 = 1-(-110 + (bias_c*5))/-90;
aj_16_p2 = 1-(-111.4 + (bias_c*5))/-90;
aj_16_p3 = 1-(-109.3 + (bias_c*5))/-90;
aj_16_p4 = 1-(-107.2 + (bias_c*5))/-90;
aj_16_p5 = 1-(-111.1 + (bias_c*5))/-90;
add_scale_n = (aj_16_n1+aj_16_n2+aj_16_n3+aj_16_n4+aj_16_n5)/5;
add_scale_p = (aj_16_p1+aj_16_p2+aj_16_p3+aj_16_p4+aj_16_p5)/5;
add_scale_final = (add_scale_n+add_scale_p)/2;
whole = [aj_16_n1 aj_16_n2 aj_16_n3 aj_16_n4 aj_16_n5
         aj_16_p1 aj_16_p2 aj_16_p3 aj_16_p4 aj_16_p5];
