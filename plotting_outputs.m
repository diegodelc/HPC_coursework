clear
close all

% read data
path = "./cw_workspace/output.txt";
data = importdata(path);

X = data(:,1);
Y = data(:,2);

u = data(:,3);
v = data(:,4);
h = data(:,5);

% plot data
figure()
surf(X,Y,h);
