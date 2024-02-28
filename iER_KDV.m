clc;
close all;
clear all;
addpath('utils\');
addpath('data\');
load SimulationData_V_1.mat;

%%
[n,m] = size(u);
[X, T] = meshgrid(x, tspan);
u = real(u);
%%
figure1 = figure
axes1 = axes('Parent',figure1);
surf(X,T,u)
map=mymap("coolwarm");
colormap(map)
%axis off
grid on
shading interp
set(axes1,'FontName','Times New Roman');
print('iER_KDV','-depsc','-vector');
%%
u = reshape(u,[n*m,1]);
ut = reshape(ut,[n*m,1]);
ut = real(ut);
X_data = [u,ut];
ux = real(ux);
uxx = real(uxx);
uxxx = real(uxxx);
ux = reshape(ux,[n*m,1]);
uxx = reshape(uxx,[n*m,1]);
uxxx = reshape(uxxx,[n*m,1]);
X_ders = [ones(n*m,1),ux,uxx,uxxx];
P = 2;
[Theta] = build_Theta(X_data, X_ders, P);

%%
% n = n/2;
% m = 2*m;
for i = 1:n
    Phi_group(:,:,i) = Theta((i-1)*m+1:i*m,:);
    ut_group(:,i) = ut((i-1)*m+1:i*m);
end

%%
clear Xi
tol = 0.01;
for j = 1:n
    for i = 1:size(Theta,2)
        B = Phi_group(:,:,j);
        B(:,[i]) = [];
        Xi(:,i,j) = erfit(B,Phi_group(:,i,j),tol);
%     err(i) = norm(Theta(:,i)-B*Xi(:,i))/norm(Theta(:,i));
    end
end
%%
clear Xi1
for k = 1:n
    Xi1(:,1,k)=[1;Xi(:,1,k)];
    for i = 1:size(Xi,2)-1
        Xi1(:,i+1,k) = [Xi(1:i,i+1,k);1;Xi(i+1:end,i+1,k)];
    end
end
%%
clear err
%uu = reshape(u,[m,n]);
for i = 1:size(Xi1,3)
    for j = 1:size(Xi1,2)
        Loss = ICcalculations(Phi_group(:,:,i),Xi1(:,j,i),ut_group(:,i));
        err(j,i)=Loss.aic_c;%norm(dxt-Theta*Xi1(:,i));

        errr(j,i) = norm(ut_group(:,i)-Phi_group(:,:,i)*Xi1(:,j,i));
    end
end
%%
[aa,bb ] = min(err);
[cc,dd ] = min(aa);
%Xi_pre = Xi1(:,bb) 
[aa1,bb1 ] = min(errr);
[cc1,dd1 ] = min(aa1);
%%
[aaa,bbb] = min(min(err))
%%
Xi_pre = Xi1(:,bb1(dd1),dd1)

%%
close all
plot(err,'ob')
hold on
plot(9,cc,'or',LineWidth=1.5)
grid on
shading interp
ylabel('AIC')
%%
figure
plot(err,'ob')
hold on
plot(9,cc,'or',LineWidth=1.5)
grid on
shading interp
ylabel('AIC','FontSize',14,'FontName','Times New Roman')
xlabel('Group','FontSize',14,'FontName','Times New Roman')
print('iERKDV_AIC','-depsc','-vector');