clc;
clear all;
close all;
addpath('.\utils');

Beta = [10;28;8/3];
x0 = [-8;8;27];
dt = 0.001;
tspan = dt:dt:50+2*dt;
[t,x] = ode45(@(t,x) Lorenz(t,x,Beta),tspan,x0);
% x0 = xt(end,:);
% [tt,x] = ode45(@(t,x) Lorenz(t,x,Beta),0:dt:1.5+2*dt,x0);
%%
plot3(x(:,1),x(:,2),x(:,3),LineWidth=1.5);
grid on
xlabel('x')
ylabel('y')
zlabel('z')
%%
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
axis off
hold(axes1,'on');

% 创建 plot3
plot3(x(:,1),x(:,2),x(:,3),'LineWidth',1.5);

view(axes1,[25.0999998771421 7.17782108672451]);
hold(axes1,'off');
%%  求导数
% 添加噪音
eps1 = 0;%0.01;
%eps2 = 0.01;
x = x + eps1*randn(size(x));

%% 使用中心差分法
dx = (1/(2*dt))*(x(3:end,:)-x(1:end-2,:));
%dx = dx + eps2*randn(size(dx));
x(1,:) = [];
x(end,:) = [];

%%  求取Phi库
% SINDy 方法确定Theta库
polyorder = 3;
%n = 3;
Theta = pooldata(0,x,polyorder,0);
%%
%组稀疏
m = 25;
n = 2000;
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
    Losses(i) = Loss.aic;
end
%%
[a,b ] = min(Losses);
Xi_pre = Xi(:,:,b)

%%

plot(Losses,'ob')
hold on
plot(b,a,'or',LineWidth=1.5)
grid on
shading interp
ylabel('AIC')




