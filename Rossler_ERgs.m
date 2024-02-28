% Copyright 2016, All Rights Reserved
% Code by Steven L. Brunton
clear all, close all, clc
figpath = './figures/';
addpath('./utils');

% generate Data
a = .2;  
b = .2;
c = 5.7;
n = 3;
x0=[1; 1; 1];  % Initial condition

% Integrate
dt = 0.01;
tspan=[dt:dt:500];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,xdat]=ode45(@(t,x) rossler(t,x,a,b,c),tspan,x0,options);

stackmax = 100;  % the number of shift-stacked rows
rmax = 6;  % maximum number of singular values
%%
figure1 = figure
axes1 = axes('Parent',figure1);
plot3(xdat(:,1),xdat(:,2),xdat(:,3))
grid on
xlabel('x','FontAngle','italic','FontSize',14,'FontName','Times New Roman')
ylabel('y','FontAngle','italic','FontSize',14,'FontName','Times New Roman')
zlabel('z','FontAngle','italic','FontSize',14,'FontName','Times New Roman')
set(axes1,'FontName','Times New Roman');
%print('rossler','-depsc','-vector');
%%
% 使用中心差分法
dx = (1/(2*dt))*(xdat(3:end,:)-xdat(1:end-2,:));
%dx = dx + eps2*randn(size(dx));
xdat(1,:) = [];
xdat(end,:) = [];
%%  求取Phi库
% SINDy 方法确定Theta库
polyorder = 2;
%n = 3;
Theta = pooldata(0,xdat,polyorder,0);
%%
%组稀疏
m = 30;
n = 1000;
clear Theta_group dx_group
for i = 1:m
    Theta_group(:,:,i) = Theta((i-1)*n+1:i*n,:);
    dx_group(:,:,i) = dx((i-1)*n+1:i*n,:);
end

%%
clear Xi
tol = 0.01;
for i = 1:m
    Xi(:,:,i) = erfit(Theta_group(:,:,i),dx_group(:,:,i),tol);
end

%%
for i = 1:m
    Loss = ICcalculations(Theta_group(:,:,i),Xi(:,:,i),dx_group(:,:,i));
    Losses(i) = Loss.aic_c;
end
%%
[a,b ] = min(Losses);
Xi_pre = Xi(:,:,b)
%%
figure
plot(Losses,'ob')
hold on
plot(b,a,'or',LineWidth=1.5)
grid on
shading interp
ylabel('AIC','FontSize',14,'FontName','Times New Roman')
xlabel('Group','FontSize',14,'FontName','Times New Roman')
print('rossler_AIC','-depsc','-vector');

