function [P00_z,K_z,z_update,zxm_z]=Kalman_Filter_update_mult_radio1070319(P00_z,Rm_data,H,R,Rm_h,k,c)

zxm_z = Rm_data(k,1:4) - Rm_h(k,1:4);

Mu_z= zxm_z';

z_update = zeros(8,1);

K_z=P00_z*H'/(H*P00_z*H'+R);    
z_update = K_z*Mu_z;
%xz_h=xz_h+z_update;                        % x(k|k)
I=eye(8);
P00_z=(I-K_z*H)*P00_z;                     % P(k|k)

%% Discrete - sequential operation
%     if k>1500
%         for i = 1:4
%             if abs(Mu_z(i,1)) < c
% %                 A=(H(i,:)*P00_z*H(i,:)'+R(i,i))
%                 K_z = P00_z*H(i,:)'/(H(i,:)*P00_z*H(i,:)'+R(i,i));
%                 I=eye(8);
%                 P00_z=(I-K_z*H(i,:))*P00_z;
%             else
%                 K_z = zeros(8,1);
%                 I=eye(8);
%                 P00_z=(I-K_z*H(i,:))*P00_z;
%             end
% 
%             z_update = z_update+K_z*Mu_z(i,1);
% 
%         end
% 
%     else
%   
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         z_update = K_z*Mu_z;
%         I=eye(8);
%         P00_z=(I-K_z*H)*P00_z;    
%     end

%% General method
% if k>500
%     
%     if (abs(zxm_z1(k)) > c && abs(zxm_z2(k)) > c && abs(zxm_z3(k)) > c && abs(zxm_z4(k)) > c) 
%         K_z = zeros(8,4);
%        Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%        Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%        Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z1(k)) > c && abs(zxm_z2(k)) > c && abs(zxm_z3(k)))  % Not Anchor1 Anchor2 Anchor3
%         Mu_z(1,:)=[];
%         Mu_z(1,:)=[];
%         Mu_z(1,:)=[];    
%         H(1,:)=[];
%         H(1,:)=[];
%         H(1,:)=[];
%         R(4,:)=[];
%         R(3,:)=[]; 
%         R(2,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         R(:,2)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=1;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z1(k)) > c && abs(zxm_z2(k)) > c && abs(zxm_z4(k))) % Not Anchor1 Anchor2 Anchor4
%         Mu_z(1,:)=[];
%         Mu_z(1,:)=[];
%         Mu_z(2,:)=[];    
%         H(1,:)=[];
%         H(1,:)=[];
%         H(2,:)=[];
%         R(4,:)=[];
%         R(3,:)=[]; 
%         R(2,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         R(:,2)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=1;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z1(k)) > c && abs(zxm_z3(k)) > c && abs(zxm_z4(k))) % Not Anchor1 Anchor3 Anchor4
%         Mu_z(1,:)=[];
%         Mu_z(2,:)=[];
%         Mu_z(2,:)=[];    
%         H(1,:)=[];
%         H(2,:)=[];
%         H(2,:)=[];
%         R(4,:)=[];
%         R(3,:)=[]; 
%         R(2,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         R(:,2)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=1;Mu_count_34=0;
%     elseif (abs(zxm_z2(k)) > c && abs(zxm_z3(k)) > c && abs(zxm_z4(k))) % Not Anchor2 Anchor3 Anchor4
%         Mu_z(2,:)=[];
%         Mu_z(2,:)=[];
%         Mu_z(2,:)=[];    
%         H(2,:)=[];
%         H(2,:)=[];
%         H(2,:)=[];
%         R(4,:)=[];
%         R(3,:)=[]; 
%         R(2,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         R(:,2)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=1;
%     elseif (abs(zxm_z1(k)) > c && abs(zxm_z2(k)) > c ) % not anchor1 anchor2
%         Mu_z(1,:)=[];
%         Mu_z(1,:)=[];
%         H(1,:)=[];
%         H(1,:)=[];
%         R(4,:)=[];
%         R(3,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_21=1;
%         i_21=i_21+1;
%         Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z1(k)) > c && abs(zxm_z3(k)) > c ) % not anchor1 anchor3
%         Mu_z(1,:)=[];
%         Mu_z(2,:)=[];
%         H(1,:)=[];
%         H(2,:)=[];
%         R(4,:)=[];
%         R(3,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_22=1;
%         i_22=i_22+1;
%         Mu_count_21=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z1(k)) > c && abs(zxm_z4(k)) > c ) % not anchor1 anchor4
%         Mu_z(1,:)=[];
%         Mu_z(3,:)=[];
%         H(1,:)=[];
%         H(3,:)=[];
%         R(4,:)=[];
%         R(3,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_23=1;
%         i_23=i_23+1;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z2(k)) > c && abs(zxm_z3(k)) > c ) % not anchor2 anchor3
%         Mu_z(2,:)=[];
%         Mu_z(2,:)=[];
%         H(2,:)=[];
%         H(2,:)=[];
%         R(4,:)=[];
%         R(3,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_24=1;
%         i_24=i_24+1;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z2(k)) > c && abs(zxm_z4(k)) > c ) % not anchor2 anchor4
%         Mu_z(2,:)=[];
%         Mu_z(3,:)=[];
%         H(2,:)=[];
%         H(3,:)=[];
%         R(4,:)=[];
%         R(3,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);    
%         Mu_count_25=1;
%         i_25=i_25+1;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_26=0;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z3(k)) > c && abs(zxm_z4(k)) > c ) % not anchor3 anchor4
%         Mu_z(3,:)=[];
%         Mu_z(3,:)=[];
%         H(3,:)=[];
%         H(3,:)=[];
%         R(4,:)=[];
%         R(3,:)=[];
%         R(:,4)=[];
%         R(:,3)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R); 
%         Mu_count_26=1;
%         i_26=i_26+1;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z1(k)) > c ) 
%         Mu_z(1,:)=[];
%         H(1,:)=[];
%         R(4,:)=[];
%         R(:,4)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_1=1;
%         i_1=i_1+1;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z2(k)) > c ) 
%         Mu_z(2,:)=[];
%         H(1,:)=[];
%         R(4,:)=[];
%         R(:,4)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R); 
%         Mu_count_2=1;
%         i_2=i_2+1;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_3=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z3(k)) > c ) 
%         Mu_z(3,:)=[];
%         H(3,:)=[];
%         R(4,:)=[];
%         R(:,4)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_3=1;
%         i_3=i_3+1;
%         Mu_count_5=0;Mu_count_4=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     elseif (abs(zxm_z4(k)) > c ) 
%         Mu_z(4,:)=[];
%         H(4,:)=[];
%         R(4,:)=[];
%         R(:,4)=[];
%         K_z = P00_z*H'/(H*P00_z*H'+R); 
%         Mu_count_4=1;
%         i_4=i_4+1;
%         Mu_count_5=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     else
%         K_z = P00_z*H'/(H*P00_z*H'+R);
%         Mu_count_5=1;
%         i_5 = i_5+1;
%         Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%         Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%         Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
%     end
% else
%     K_z = P00_z*H'/(H*P00_z*H'+R);    
%     Mu_count_5=1;
%     i_5 = i_5+1;
%     Mu_count_4=0;Mu_count_3=0;Mu_count_2=0;Mu_count_1=0;
%     Mu_count_21=0;Mu_count_22=0;Mu_count_23=0;Mu_count_24=0;Mu_count_25=0;Mu_count_26=0;
%     Mu_count_31=0;Mu_count_32=0;Mu_count_33=0;Mu_count_34=0;
% 
% end



% %%
% 
% z_update = K_z*Mu_z;
% I=eye(8);
% P00_z=(I-K_z*H)*P00_z;                     % P(k|k)
% 
% end

end

