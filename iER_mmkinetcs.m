clear all, close all, clc
addpath('./utils');

%% Define System

%size of system
n = 1;
% rate constants
Vmax = 1.5; % maximum reaction rate
Km = 0.3;   % half max reaction rate
jin = 0.6;  % influx of substrate

% reaction function
MMKinetics = @(x)jin-Vmax.*x./(Km+x);
%% generate Data
measure = 2; % number of initial conditions
dt = 0.1; % time step saved
tspan=[0:dt:4]; % time vector
N = length(tspan);

% intial condition of concentration
x0 = 0.5;

for ii = 1:measure
    % Integrate
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    [t1,x1]=ode45(@(t,x)MMKinetics(x),tspan,x0,options);
    tt(:,ii) = t1;
    x(:,:,ii) = x1;
    x0 = 2*x0; % get a new initial condition
end

%% add noise & calculate derivative
eps = 1e-4; %magnitude of noise
xn = x + eps*randn(size(x)); % add normally distributed measuement error

%calculate exactly and add error
xt = [];dxt= []; t = [];

for ll =1:measure
    for ii=1:length(x)
        dxf(ii,:,ll) = MMKinetics(x(ii,:, ll));
    end
    epsdt =0;
    dxf = dxf+epsdt*randn(size(dxf));
    dxt = [dxt; dxf(:,:,ll)];
    xt = [xt; xn(:,:,ll)];
    t = [t; tt(:, ll)];
end
% figure(5)
% hold off
% plot(t ,xt, 'o')
% xlabel('time')
% ylabel('concentrations')
% title('training time series')
% 
% figure(6) 
% hold off
% plot(t, dxt, 'o')
% xlabel('time')
% ylabel('derivative of concentrations w/ time')
% title('training derivative time series')

%% set the functions to search
% may not be the same as the functions we used to generate the data
laurentorder = 0; % n, for 1/x^n terms in the library
polyorder = 4; % n, for polynomial terms in the library x^n
usesine = 0; % adds sine to the library
dyorder = 1; % n for (dx/dt)^n terms in the library

% pool Data  (i.e., build library of nonlinear time series)
[Theta, Thetastring] = poolDatady(xt,n,polyorder,usesine, laurentorder, dxt, dyorder);

%%
tol = 0.001;
for i = 1:size(Theta,2)
    B = Theta;
    B(:,[i]) = [];
    Xi(:,i) = erfit(B,Theta(:,i),tol);
%     err(i) = norm(Theta(:,i)-B*Xi(:,i))/norm(Theta(:,i));
end
%%
clear Xi1
Xi1(:,1)=[1;Xi(:,1)];
for i = 1:size(Xi,2)-1
    Xi1(:,i+1) = [Xi(1:i,i+1);1;Xi(i+1:end,i+1)];
end
%%
clear err
for i = 1:size(Xi1,2)
    Loss = ICcalculations(Theta,Xi1(:,i),dxt);
    err(i)=Loss.aic_c;%norm(dxt-Theta*Xi1(:,i));
end

%%
[aa,bb ] = min(err)
Xi_pre = Xi1(:,bb) 
%%
figure1 = figure
axes1 = axes('Parent',figure1);
hold off
plot(t ,xt, 'ob')
hold on
plot(t, dxt, 'or')
grid on
legend('x','derivative','FontSize',14,'FontName','Times New Roman')
set(axes1,'FontName','Times New Roman');
print('MM','-depsc','-vector');
%%
figure
plot(err,'ob')
hold on 
plot(bb,aa,'or',LineWidth=2)
ylabel('AIC','FontSize',14,'FontName','Times New Roman')
xlabel('Group','FontSize',14,'FontName','Times New Roman')
grid on
print('MMAIC','-depsc','-vector');


