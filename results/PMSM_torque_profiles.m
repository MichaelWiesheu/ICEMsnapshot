clc
clear all
close all

restoredefaultpath
addpath(genpath("../"))

%%
figure(1)
clf
load('InitTorques.mat');
angles = 0:20;
Torques = [Torques, Torques(:,1)];
lw = 1.5;
hold on
plot(angles, Torques(1,:), "Color", TUDa_getColor_num(1), LineWidth=lw);
plot(angles, Torques(2,:), "Color", TUDa_getColor_num(1), LineWidth=lw);
plot(angles, Torques(3,:), "Color", TUDa_getColor_num(1), LineWidth=lw);

data(:,1) = angles;
dataScatter(:,1) = angles(1:2:end-1);
data(:,2:4) = Torques';
dataScatter(:,2:4) = Torques(:,1:2:end-1)';
%%
load('OptTorques.mat')
angles = 0:20;
Torques = [Torques, Torques(:,1)];

hold on
plot(angles, Torques(1,:), "Color", TUDa_getColor_num(3), LineWidth=lw);
plot(angles, Torques(2,:), "Color", TUDa_getColor_num(3), LineWidth=lw);
plot(angles, Torques(3,:), "Color", TUDa_getColor_num(3), LineWidth=lw);

grid on

xlabel("Rotation angle (°)")
ylabel("Torque (Nm)")

data(:,5:7) = Torques';
dataScatter(:,5:7) = Torques(:,1:2:end-1)';

%%

matlab2tikz('figurehandle', gcf, 'Torques.tikz')

csvwrite('TorqueData.csv', data)
csvwrite('TorqueScatter.csv', dataScatter)


