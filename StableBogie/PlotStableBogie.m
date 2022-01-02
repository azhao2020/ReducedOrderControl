close all
clear all


swEPSfigure

load StableBogie.mat
%%

%=========================================================%
% Plot the results from Simulink
% Figure(1): Control and controlled system state
% Figure(2): Unforced original and reduced system
%=========================================================%
figure
swFigSize
subplot(2,1,1)
plot(tspan(100:100000), x_o_u0(100:100000,1),'--r','LineWidth',1)
hold on
plot(tspan(1:100000), x_o_u(1:100000,1),'-k','LineWidth',1)
% xlabel('Time ($s$)')
ylabel('$x_1(t) (m)$')
axis([0 6 -0.05 0.05])
legend('Open-loop','Closed-loop')

subplot(2,1,2)
plot(tspan(100:30000),u(100:30000,1),'--r','LineWidth',1);
hold on
plot(tspan(100:30000),u(100:30000,2),'-k','LineWidth',1);
xlabel('Time ($s$)')
ylabel('Control (N)')
legend('$u_1(t)$','$u_2(t)$','Interpreter','latex')
% axis([0.1 3 -20 20])
print -depsc SystemStableControlled.eps

% figure(2)
% i = length(tspan);
% for ite = 1:i
%     norm_yr(ite) = norm(y_o_u0(ite,:));
%     norm_yo(ite) = norm(y_r_u0(ite,:));
% end
% plot(tspan,norm_yo,'--');
% hold on
% plot(tspan,norm_yr)
% swFigSize
% xlabel('Time ($s$)')
% ylabel('Output Norm')
% axis([0 10 0 0.06])
% legend('$||\mathbf{y}_o(t)||$','$||\mathbf{y}_r(t)||$','Interpreter','latex')
% print -depsc SystemOutputFree.eps
% 
% %================================================%
% % Plot norm of the error between the outputs
% %================================================%
% figure(3)
% plot(tspan,norm_yo-norm_yr);
% swFigSize
% xlabel('Time ($s$)')
% ylabel('Output Error Norm')
% axis([0 10 0 0.001])
% legend('$||\mathbf{e}_{y(t)}||$','Interpreter','latex')
% print -depsc SystemOutputError.eps



figure(2)
subplot(2,1,1)
i = length(tspan);
for ite = 1:i
    norm_yr(ite) = norm(y_o_u0(ite,:));
    norm_yo(ite) = norm(y_r_u0(ite,:));
end
plot(tspan,norm_yo,'--');
hold on
plot(tspan,norm_yr)
% swFigSize
% xlabel('Time ($s$)')
ylabel('Output Norm')
axis([0 4 0 0.05])
legend('$||\mathbf{y}_o(t)||$','$||\mathbf{y}_r(t)||$','Interpreter','latex')



%================================================%
% Plot norm of the error between the outputs
%================================================%
% figure(3)
subplot(2,1,2)
norm_error = norm_yo-norm_yr;
plot(tspan,norm_error);
% swFigSize
xlabel('Time ($s$)')
ylabel('$||\mathbf{e}_{y(t)}||_{L_2}$')
axis([0 4 0 0.001])
legend('$||\mathbf{e}_{y(t)}||$','Interpreter','latex')
print -depsc SystemOutputFree.eps


%================================================%
% Plot Hankel singular values for balanced system
% Different illustration figures
%================================================%
hsvs = hsvd(sysb);
% figure()
% set(gcf,'Position',[100 100 800 700]);
% h1 = hsvplot(sysb); 
% print -depsc HSVStable.eps

% figure()
% plot(hsvs,'o-')

figure(4)
bar_h1 = semilogy([1:length(hsvs)], hsvs, 'o-b');
swFigSize
set(gca,'XTick',[1:2:16])
xlabel('State')
ylabel('Hankel Singular Value $\sigma_i$')
legend('Stable modes')
print -depsc HSVStable.eps

%=============================================================%
% Numerical criteria to describe the perfomance of the control
%=============================================================%
% Criterion1: Using ||x||
for i = 1:size(x_o_u,1)
    x_o_norm(i) = norm(x_o_u(i,:));
    x_o_u0_norm(i) = norm(x_o_u0(i,:));
end
    
index1 = find(x_o_norm<0.01);
time1 = index1(1)*10^-4;
text1 = sprintf('It takes %d seconds for the closed-loop system in stable case to reach ||x||<0.01',time1);
disp(text1)

index2 = find(x_o_u0_norm<0.01);
time2 = index2(1)*10^-4;
text2 = sprintf('It takes %d seconds for the open-loop system in stable case to reach ||x||<0.01',time2);
disp(text2)

figure(5)
plot(tspan,x_o_u0_norm,'-.r','LineWidth',0.5)
hold on
plot(tspan,x_o_norm,'-k','LineWidth',1)
axis([0 6 -0.1 0.1])

xlabel('Time ($s$)')
ylabel('$||\mathbf{x}(t)||_{L_2}$')
legend('Open-loop','Closed-loop')
a=axis;
axis([0, 6, 0, 1])
print -depsc StateNormStableCase.eps

%===================================%
% Estimated state and error
%===================================%
subplot(2,1,1)
plot(tspan(100:30000), x_r_u(100:30000,1))
xlabel('Time ($s$)')
ylabel('$\hat{x}_{r1}(t)$')

subplot(2,1,2)
for i = 100:30000
     norm_yest_err(i) = norm(y_est_err(i,:));
end
plot(tspan(100:30000),norm_yest_err(100:30000))
xlabel('Time ($s$)')
ylabel('$||\mathbf{e}_{y}(t)||$')

print -depsc EstimatedStateandError