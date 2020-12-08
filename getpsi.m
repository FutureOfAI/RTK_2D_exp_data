close all; %close all figures
clear all;
clc;       %clear the command line
fclose('all'); %close all open files
delete(instrfindall); %Reset Com Port

obj = serial('com10');
obj.BaudRate = 115200;
obj.Parity = 'none';
obj.StopBits = 1;
count = 2000;
fopen(obj);
disp(obj);
data = zeros(5,count);
for i=1:count
    data(:,i) = fread(obj,5);
end

fclose(obj);
delete(obj);

save aj_psi.mat data;