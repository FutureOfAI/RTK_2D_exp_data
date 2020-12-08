%% data format
% ax ay az in g
% gx gy gz in deg/s
% longitude latitude (in deg) height(in m)
% vel_east vel_north in m/s
% odometer left/right wheel in m
% V_odo left/right wheel in m/s
% dt in s

format long

% data = cell2mat(struct2cell(load('2019-10-15-14-47-28am.mat'))); % load data to matrix
data = cell2mat(struct2cell(load('2019-10-15-16-37-46pm.mat'))); % load data to matrix
attu = cell2mat(struct2cell(load('robotAttitude.mat')));

[row, col] = size(data);
t = 1: row;

%% BLH-frame to ECEF-frame
earthR = 6378137.0; % earth radius in m
earth_W = 0.000072921151467; % rad/s
earth_F = 1/298.257223563;
earth_E2 = earth_F*(2-earth_F);
r2d = (180/pi);
d2r = (pi/180);
g = 9.8;
acc_g = [0 0 g];
corECEF = zeros(row,3);
p_N = zeros(row,3);
e2 = earth_F*(2-earth_F); % ÍÖÇò±âÂÊ
station_p = [-2810989.83000001 4717661.46600002 3233015.55499995]; % station position

for i = 1:row
    phi = data(i, 7)*d2r; % lat in rad
    lambda = data(i, 8)*d2r; % long in rad
    h = data(i, 9);
    N = earthR/sqrt(1-e2*sin(phi)*sin(phi));
    corECEF(i,:) = [(N+h)*cos(phi)*cos(lambda) (N+h)*cos(phi)*sin(lambda) (N*(1-e2)+h)*sin(phi)];
end

% defferential GPS
detaP = zeros(row,3);
for i = 1:row
    detaP(i, :) = corECEF(i, :) - station_p;
end

% ECF to E-N-U
conversionMatrix = zeros(3,3);
for i = 1:row
    lat = data(i,7)*d2r;
    long = data(i,8)*d2r;
    C_ECF2N = cz(pi/2)*cy(-lat+pi/2)*cz(long);
    p_N(i,:) = (C_ECF2N*detaP(i,:)')'; % antana
end

attz = zeros(row, 1);
for i = 1:row
    if attu(i,3)>pi/2
        attz(i) = attu(i,3)-2*pi;
    else
        attz(i) = attu(i,3);
    end
end

accx = data(t,1)*g;
accy = data(t,2)*g;
accz = data(t,3)*g;
grox = data(t,4)*d2r;
groy = data(t,5)*d2r;
groz = data(t,6)*d2r;
