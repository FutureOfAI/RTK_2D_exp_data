close all; %close all figures
clear all;
clc;       %clear the command line
fclose('all'); %close all open files
delete(instrfindall); %Reset Com Port

obj = serial('com10');
obj.BaudRate = 115200;
obj.Parity = 'none';
obj.StopBits = 1;
count = 200;
fopen(obj);
disp(obj);
data = zeros(9,count);
for i=1:count
    data(:,i) = fread(obj,9);
end

fclose(obj);
delete(obj);

save position1.mat data;