clc;
close all;
clear all;
addpath('utils\');
addpath('data\');
load burgers.mat;
U = real(usol);
[n,m] = size(U);
[X, T] = meshgrid(x, t);
%U = U+0.01*std2(U)*randn(size(U));
%%
figure1 = figure
axes1 = axes('Parent',figure1);
surf(X,T,U')
map=mymap("coolwarm");
colormap(map)
%axis off
grid on
shading interp
set(axes1,'FontName','Times New Roman');
print('Burger','-depsc','-vector');
%%
dt = t(2)-t(1);
dx = x(2)-x(1);
D = 3;
P = 3;

for i = 1:n
    ut(i,:) = derivative(U(i,:),dt,2);
end
%%
ut = reshape(ut,[n*m,1]);

%% 构建Phi库
Phi = ones(n*m,(D+1)*(P+1));
ux = zeros(n,m);
for d = 1:D+1
    if d>1
        for i = 1:m
            ux(:,i) = derivative(U(:,i),dx,d);
        end
    else
        ux = ones(n,m);
    end 
        for p = 1:P+1
            u1 = ux.*(U.^(p-1));
            Phi(:,(d-1)*(P+1)+p) = reshape(u1,[n*m,1]);
        end

end
%%
%组稀疏
m = 101;
n = 236;
clear Phi_group ut_group
for i = 1:m
    Phi_group(:,:,i) = Phi((i-1)*n+1:i*n,:);
    ut_group(:,i) = ut((i-1)*n+1:i*n);
end
%%
tol = 0.05;
for i = 1:m
    Xi(:,i) = erfit(Phi_group(:,:,i),ut_group(:,i),tol);
end
%%
for i = 1:m
    Loss = ICcalculations(Phi_group(:,:,i),Xi(:,i),ut_group(:,i));
    Losses(i) = Loss.aic;
end
%%
[a,b ] = min(Losses);
Xi_pre = Xi(:,b)
% lam = 1e-5;
% d_tol = 3;
% [w_best,tol_best] = TrainSTRidge(Phi, ut, lam, d_tol)
%%
figure
plot(Losses,'ob')
hold on
plot(b,a,'or',LineWidth=1.5)
grid on
shading interp
ylabel('AIC','FontSize',14,'FontName','Times New Roman')
xlabel('Group','FontSize',14,'FontName','Times New Roman')
print('Burger_AIC','-depsc','-vector');
