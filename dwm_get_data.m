close all; %close all figures
clear all;
clc;       %clear the command line
fclose('all'); %close all open files
delete(instrfindall); %Reset Com Port

obj = serial('com18');
obj.BaudRate = 57600;
obj.Parity = 'none';
obj.StopBits = 1;
count = 4000;
fopen(obj);
disp(obj);
data = zeros(44,count);
for i=1:count
    data(:,i) = fread(obj,44);
end

fclose(obj);
delete(obj);

save multisapce_data_1070305_v2.mat data;