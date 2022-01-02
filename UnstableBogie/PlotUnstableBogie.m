close all
clear all

swEPSfigure

load UnstableBogie.mat
%%
%=========================================================%
% Plot the results from Simulink
% Figure(1): Control Torque
% Figure(2): Unforced original and reduced system
%=========================================================%
% swFigSize
figure(1)
% swFigSize
subplot(2,1,1)
plot(tspan(100:80000), x_o_u(100:80000,1),'-k','LineWidth',1)
% xlabel('Time ($s$)')
ylabel('$x_1(t) (m)$')
% legend('$x_1$','Interpreter','latex')

subplot(2,1,2)
plot(tspan(100:20000),u(100:20000,1),'--r','LineWidth',1);
hold on
plot(tspan(100:20000),u(100:20000,2),'-k','LineWidth',1);
xlabel('Time ($s$)')
ylabel('Control (N)')
legend('$u_1(t)$','$u_2(t)$','Interpreter','latex')

print -depsc SystemUnstableControlled.eps

figure(2)
subplot(2,1,1)
i = length(tspan);
for ite = 1:i
    norm_yr(ite) = norm(y_o_u0(ite,:));
    norm_yo(ite) = norm(y_r_u0(ite,:));
end
plot(tspan(1:40000,:),norm_yo(:,1:40000),'--');
hold on
plot(tspan(1:40000,:),norm_yr(:,1:40000))
% swFigSize
% xlabel('Time ($s$)')
ylabel('Output Norm')
axis([0 4 0 0.15])
legend('$||\mathbf{y}_o(t)||$','$||\mathbf{y}_r(t)||$','Interpreter','latex')

%================================================%
% Plot norm of the error between the outputs
%================================================%
% figure(3)
subplot(2,1,2)
norm_error = norm_yo-norm_yr;
plot(tspan(1:40000,:),norm_error(:,1:40000));
% swFigSize
xlabel('Time ($s$)')
ylabel('$||\mathbf{e}_{y(t)}||_{L_2}$')
axis([0 4 0 0.001])
legend('$||\mathbf{e}_{y(t)}||$','Interpreter','latex')
print -depsc SystemUnstableOutputFree.eps

%================================================%
% Plot Hankel singular values for balanced system
%================================================%
hsvs = hsvd(sysb);
idx = find(hsvs==Inf);
hsvs(idx) = 0.25;
hsvs1 = hsvs(hsvs>0.2);
hsvs2 = hsvs(hsvs<0.2);

% bar plot
figure(4)
swFigSize
bar_h1 = semilogy([1:length(idx)], hsvs1,'-*r',[length(idx)+1:length(hsvs)], hsvs2, 'o-b');
set(gca,'XTick',[1:2:16])
xlabel('State')
ylabel('Hankel Singular Value $\sigma_i$')
legend('Unstable states','Stable states', 'Location','NE')
print -depsc HSVUnstable.eps

%=============================================================%
% Numerical criteria to describe the perfomance of the control
%=============================================================%
% Criterion1: Using ||x||
for i = 1:size(x_o_u,1)
    x_o_norm(i) = norm(x_o_u(i,:));
end
    
index1 = find(x_o_norm<0.01);
time1 = index1(1)*10^-4;
text1 = sprintf('It takes %d seconds for the stable system to reach ||x||<0.01',time1);
disp(text1)

figure(5)
swFigSize
plot(tspan,x_o_norm,'-k','LineWidth',1)

xlabel('Time ($s$)')
ylabel('$||\mathbf{x}(t)||_{L_2}$')
axis([0 8 0 1.3])
legend('Closed-loop')

print -depsc StateNormUnstableCase.eps

%%
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

print -depsc EstimatedStateandErrorUnstableCase.eps