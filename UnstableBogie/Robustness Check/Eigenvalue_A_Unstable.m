%=============================================================%
% Eigenvalue changing of the system
% Robustness check
%=============================================================%

clc
clear
close all
format long
% swEPSfigure

load UnstableBogie.mat
%%
ite = 0;
for v = 50:170
    ite = ite + 1;
    w = 1/v;
    CD = [2*fy*w,0,0,0,0,0,0,0;
        0,2*l0^2*fx*w,0,0,0,0,0,0;
        0,0,2*fy*w,0,0,0,0,0;
        0,0,0,2*l0^2*fx*w,0,0,0,0;
        zeros(4,8)];
    A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
    
%   eig_Acl(ite) = max(real(eig(A-Bo*(K1*inv(C_r)*C))));
%   eig_Acl(ite) = max(real(eig(A-Bo*(K1*Tr))));
    
    zeta0(ite) = min(-real(eig(A-Bo*(K1*Tr)))/(norm(eig(A-Bo*(K1*Tr)))));

    v_var(ite) = v*3.6;
end

idx = find(zeta0>0);
zeta0 = zeta0(idx);

swFigSize
swEPSfigure
figure(1)
plot(v_var(idx), zeta0,'-k','LineWidth',1)
xlabel('Velocity (km/h)')
ylabel('Damping ratio $\zeta$')
% %%
% % %=============================================================%
% % % Check the transformation from xo to xr 
% % % xr = T*xo
% % %=============================================================%
% % clc
% % clear
% % load UnstableBogie.mat
% 
% tspan = simOut.get('tout');
% x_o_u = simOut.get('x_ori');
% x_r_u = simOut.get('x_red');
% x_o_r = (T*x_o_u')';
% 
% figure(2)
% plot(tspan(1:200:end),x_o_r(1:200:end,1),'-k','LineWidth',1.5)
% hold on
% plot(tspan(1:200:end),x_r_u(1:200:end,1),'or')
% axis([0 8 -0.02 0.02])
% legend('Reduced states after transformation','Reduced states from simulation')
%%
v = 100;                        % Velocity in m/s
w = 1/v;
CA = [zeros(4,8);
    C2];
CD = [2*fy*w,0,0,0,0,0,0,0;
    0,2*l0^2*fx*w,0,0,0,0,0,0;
    0,0,2*fy*w,0,0,0,0,0;
    0,0,0,2*l0^2*fx*w,0,0,0,0;
    zeros(4,8)];
le = 1;
E = [0 0; 0 0; 0 0; 0 0; 1 1; le -le; 0 0; 0 0];
W = [kpy;0;kpy;0;0;0;0;0];
A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];

% Create new matrix B
B = [zeros(8,2);M\E];
Bo = B;
% z0 = norm(x_o)/norm(B(:,1));
z0 = a;
B = [B,x_o/z0];                                        % x_o/z0 is a unit vector

% Controllability check
p = 4;
q = 12;
C = [eye(p),zeros(p,q)];
D = zeros(p,size(B,2));
sys = ss(A,B,C,D);

[sysb,g,T,Ti] = balreal(sys);

order_s = 5;                                            % Truncated order
order_e = 16;
% sysmdc = modred(sysb,order_s:order_e,'MathDC');       % Reduced order subsystem
sysdel = modred(sysb,order_s:order_e,'Truncate');       % Reduced order subsystem
order = order_s-1;

Ti = Ti(:,1:order);
Tr = T(1:order,:);

% Reduced system
A_r = sysdel.A;
B_r = sysdel.B;
B_r = B_r(:,1:2);
C_r = sysdel.C;
x_r = Tr*x_o;

% Feedback gain
Q = 10^5*eye(order);
R = [10,0;0,10]/10^7;
[K1,S1,P1] = lqr(A_r,B_r,Q,R);

ite = 0;
for v = 50:170
    ite = ite + 1;
    w = 1/v;
    CD = [2*fy*w,0,0,0,0,0,0,0;
        0,2*l0^2*fx*w,0,0,0,0,0,0;
        0,0,2*fy*w,0,0,0,0,0;
        0,0,0,2*l0^2*fx*w,0,0,0,0;
        zeros(4,8)];
    A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
    
%     eig_Acl(ite) = max(real(eig(A-Bo*(K1*inv(C_r)*C))));
%    eig_Acl(ite) = max(real(eig(A-Bo*(K1*T))));

    zeta1(ite) = min(-real(eig(A-Bo*(K1*Tr)))/(norm(eig(A-Bo*(K1*Tr)))));

    v_var(ite) = v*3.6;
end


idx = find(zeta1>0);
zeta1 = zeta1(idx);

hold on
plot(v_var(idx), zeta1,'-.k','LineWidth',1.5)
xlabel('Velocity (km/h)')
ylabel('Damping ratio $\zeta$')

% Feedback gain
Q = 10^5*eye(order);
R = [10,0;0,10]/10^8;
[K1,S1,P1] = lqr(A_r,B_r,Q,R);

ite = 0;
for v = 50:200
    ite = ite + 1;
    w = 1/v;
    CD = [2*fy*w,0,0,0,0,0,0,0;
        0,2*l0^2*fx*w,0,0,0,0,0,0;
        0,0,2*fy*w,0,0,0,0,0;
        0,0,0,2*l0^2*fx*w,0,0,0,0;
        zeros(4,8)];
    A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
    
%   eig_Acl(ite) = max(real(eig(A-Bo*(K1*inv(C_r)*C))));
%   eig_Acl(ite) = max(real(eig(A-Bo*(K1*T))));

    zeta2(ite) = min(-real(eig(A-Bo*(K1*Tr)))/(norm(eig(A-Bo*(K1*Tr)))));

    v_var(ite) = v*3.6;
end

idx = find(zeta2>0);
zeta2 = zeta2(idx);

hold on
plot(v_var(idx), zeta2,'--k','LineWidth',1.5)
xlabel('Train Speed ($km/h$)')
ylabel('Damping ratio $\zeta$')



% Critical speed
for v = 50:200
    ite = ite + 1;
    w = 1/v;
    CD = [2*fy*w,0,0,0,0,0,0,0;
        0,2*l0^2*fx*w,0,0,0,0,0,0;
        0,0,2*fy*w,0,0,0,0,0;
        0,0,0,2*l0^2*fx*w,0,0,0,0;
        zeros(4,8)];
    A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
    
%   eig_Acl(ite) = max(real(eig(A-Bo*(K1*inv(C_r)*C))));
%   eig_Acl(ite) = max(real(eig(A-Bo*(K1*T))));

    zeta3(ite) = min(-real(eig(A))/(norm(eig(A))));

    v_var(ite) = v*3.6;
end

idx = find(zeta3>0);
zeta3 = zeta3(idx);

hold on
plot(v_var(idx), zeta3,':k','LineWidth',1.5)
xlabel('Train Speed ($km/h$)')
ylabel('Damping ratio $\zeta_{min}$')


leg = legend('$v_{target}=432 km/h$','$v_{target}=360 km/h$','$v_{target}=400 km/h$','open-loop','Interpreter','latex');
set(leg,'location','best')
hold on
ylim([0 7*10^-3])
xlim([100 725])
text(360, 2.8*10^-3,'$\diamondsuit$','color','b')
text(432, 1.65*10^-3,'$\diamondsuit$','color','b')
text(360, 5.1*10^-3,'$\diamondsuit$','color','b')
text(400, 5.4*10^-3,'$\nabla$','color','g')
text(415,0,'o','color','r')
text(449,0,'o','color','r')
text(491,0,'o','color','r')
text(678,0,'o','color','r')



% annotation('arrow',[0.45 0.55],[0.7 0.7],'Color','r');
% text(210,0.3,'Critical speed','FontSize',12);
% text(150,-0.2,'Stability boundary','FontSize',12);
print -depsc EigAChanging.eps
% set(gcf, 'Position', [4.874200000000000e+03 2.542000000000000e+02 560 420])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Unstable case
% lambdae-zeta relationship
% clc
% clear
% close all
% format long
% swEPSfigure
% 
% W= 21000*9.81;
% mm=3530/2;     % 1765kg                  
% Im=1019/2;
% 
% mw=2585;
% Iw=2024;
% lm=0.5;
% 
% mf=3656+mm*2;               
% If=7849+2*mm*lm^2+Im*2;        
%     
% csx=1000000*1;
% csy=60000*1;
% kpx=4*10^7;  % 4*10^7;
% l0=0.7465;
% kpy=6*10^6; %+lambda*W/l0;
% 
% kgfai=-23000;       
% ksx=0.2*10^6*4 ;   
% ksy=0.2*10^6*4 ;    
% fy=8.62400*10^6;
% fx=8.14400*10^6;
% 
% l1=1.1;
% l2=0.95;
% b=2.5/2;
% r0=0.625;
% 
% freqmy=1.8;
% dumpmy=0.2;
% 
% kym=(2*pi*freqmy)^2*mm;
% cym=dumpmy*2*mm*(2*pi*freqmy);
% M0=[mw,Iw,mw,Iw,mf,If,mm,mm];
% M=diag(M0);
% 
% lambda=0.1;
% kgy=lambda*W/l0;
% Ks = [ kpy+kgy,0,0,0,-kpy,-b*kpy,0,0;
%      0,l1^2*kpx+kgfai,0,0,0,-l1^2*kpx,0,0;
%      0,0,kpy+kgy,0,-kpy,b*kpy,0,0;
%      0,0,0,l1^2*kpx+kgfai,0,-l1^2*kpx,0,0;
%      -kpy,0,-kpy,0,(2*kpy+ksy+2*kym),0,-kym,-kym;
%      -b*kpy,-l1^2*kpx,b*kpy,-l1^2*kpx,0,l2^2*ksx+2*l1^2*kpx+2*b^2*kpy+2*lm^2*kym,-lm*kym,lm*kym;
%      0,0,0,0,-kym,-lm*kym,kym,0;   
%      0,0,0,0,-kym,lm*kym,0,kym ];
%  
% Kus = [ 0,-2*fy,0,0,0,0,0,0;
%      2*lambda*l0*fx/r0,0,0,0,0,0,0,0;
%      0,0,0,-2*fy,0,0,0,0;
%      0,0,2*lambda*l0*fx/r0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0;   
%      0,0,0,0,0,0,0,0]; 
%  
% C2 = [0,0,0,0,csy+2*cym,0,-cym,-cym;
%     0,0,0,0,0,l2^2*csx+lm^2*cym,-lm*cym,lm*cym;
%     0,0,0,0,-cym,-lm*cym,cym,0;
%     0,0,0,0,-cym,lm*cym,0,cym ];
% v = 120;                        % Velocity in m/s
% w = 1/v;
% CA = [zeros(4,8);
%     C2];
% CD = [2*fy*w,0,0,0,0,0,0,0;
%     0,2*l0^2*fx*w,0,0,0,0,0,0;
%     0,0,2*fy*w,0,0,0,0,0;
%     0,0,0,2*l0^2*fx*w,0,0,0,0;
%     zeros(4,8)];
% le = 1;
% E = [0 0; 0 0; 0 0; 0 0; 1 1; le -le; 0 0; 0 0];
% A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
% 
% With initial condition
% Set initial condition x_o(0)
% a = 0.05;                                             % Initial condition 
% x_o = a*ones(16,1);                      
% 
% Create new matrix B
% B = [zeros(8,2);M\E];
% Bo = B;
% z0 = norm(x_o)/norm(B(:,1));
% z0 = a;
% B = [B,x_o/z0];                                        % x_o/z0 is a unit vector
% 
% Controllability check
% p = 4;
% q = 12;
% C = [eye(p),zeros(p,q)];
% D = zeros(p,size(B,2));
% sys = ss(A,B,C,D);
% 
% [sysb,g,T,Ti] = balreal(sys);
% 
% order_s = 5;                                            % Truncated order
% order_e = 16;
% sysmdc = modred(sysb,order_s:order_e,'MathDC');       % Reduced order subsystem
% sysdel = modred(sysb,order_s:order_e,'Truncate');       % Reduced order subsystem
% order = order_s-1;
% 
% Ti = Ti(1:order,:);
% T = T(:,1:order);
% 
% Reduced system
% A_r = sysdel.A;
% B_r = sysdel.B;
% B_r = B_r(:,1:2);
% C_r = sysdel.C;
% x_r = Ti*x_o;
% 
% Feedback gain
% Q = 10*eye(order);
% R = [10,0;0,10]/10^7;
% [K1,S1,P1] = lqr(A_r,B_r,Q,R);
% 
% ite = 0;
% for lambda = 0.1:0.01:0.2
%     ite = ite + 1;
%     kgy=lambda*W/l0;
%     Ks = [ kpy+kgy,0,0,0,-kpy,-b*kpy,0,0;
%            0,l1^2*kpx+kgfai,0,0,0,-l1^2*kpx,0,0;
%            0,0,kpy+kgy,0,-kpy,b*kpy,0,0;
%            0,0,0,l1^2*kpx+kgfai,0,-l1^2*kpx,0,0;
%            -kpy,0,-kpy,0,(2*kpy+ksy+2*kym),0,-kym,-kym;
%            -b*kpy,-l1^2*kpx,b*kpy,-l1^2*kpx,0,l2^2*ksx+2*l1^2*kpx+2*b^2*kpy+2*lm^2*kym,-lm*kym,lm*kym;
%            0,0,0,0,-kym,-lm*kym,kym,0;   
%            0,0,0,0,-kym,lm*kym,0,kym ];
%  
%     Kus = [ 0,-2*fy,0,0,0,0,0,0;
%         2*lambda*l0*fx/r0,0,0,0,0,0,0,0;
%         0,0,0,-2*fy,0,0,0,0;
%         0,0,2*lambda*l0*fx/r0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0;   
%         0,0,0,0,0,0,0,0]; 
%     A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
%     
%     eig_A(ite) = max(real(eig(A-Bo*(K1*pinv(C_r)*C))));
%     lambda_var(ite) = lambda;
% end
% 
% swFigSize
% swEPSfigure
% plot(lambda_var, eig_A)
% xlabel('$\lambda_e$')
% ylabel('$\zeta_{max}$')
% title('Unstable case')
% set(gcf, 'Position', [4.874200000000000e+03 2.542000000000000e+02 560 420])
% %% Stable case
% % lambdae-zeta relationship
% clc
% clear
% close all
% format long
% % swEPSfigure
% 
% W= 21000*9.81;
% mm=3530/2;     % 1765kg                  
% Im=1019/2;
% 
% mw=2585;
% Iw=2024;
% lm=0.5;
% 
% mf=3656+mm*2;               
% If=7849+2*mm*lm^2+Im*2;        
%     
% csx=1000000*1;
% csy=60000*1;
% kpx=4*10^7;  % 4*10^7;
% l0=0.7465;
% kpy=6*10^6; %+lambda*W/l0;
% 
% kgfai=-23000;       
% ksx=0.2*10^6*4 ;   
% ksy=0.2*10^6*4 ;    
% fy=8.62400*10^6;
% fx=8.14400*10^6;
% 
% l1=1.1;
% l2=0.95;
% b=2.5/2;
% r0=0.625;
% 
% freqmy=1.8;
% dumpmy=0.2;
% 
% kym=(2*pi*freqmy)^2*mm;
% cym=dumpmy*2*mm*(2*pi*freqmy);
% M0=[mw,Iw,mw,Iw,mf,If,mm,mm];
% M=diag(M0);
% 
% lambda=0.1;
% kgy=lambda*W/l0;
% Ks = [ kpy+kgy,0,0,0,-kpy,-b*kpy,0,0;
%      0,l1^2*kpx+kgfai,0,0,0,-l1^2*kpx,0,0;
%      0,0,kpy+kgy,0,-kpy,b*kpy,0,0;
%      0,0,0,l1^2*kpx+kgfai,0,-l1^2*kpx,0,0;
%      -kpy,0,-kpy,0,(2*kpy+ksy+2*kym),0,-kym,-kym;
%      -b*kpy,-l1^2*kpx,b*kpy,-l1^2*kpx,0,l2^2*ksx+2*l1^2*kpx+2*b^2*kpy+2*lm^2*kym,-lm*kym,lm*kym;
%      0,0,0,0,-kym,-lm*kym,kym,0;   
%      0,0,0,0,-kym,lm*kym,0,kym ];
%  
% Kus = [ 0,-2*fy,0,0,0,0,0,0;
%      2*lambda*l0*fx/r0,0,0,0,0,0,0,0;
%      0,0,0,-2*fy,0,0,0,0;
%      0,0,2*lambda*l0*fx/r0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0;   
%      0,0,0,0,0,0,0,0]; 
%  
% C2 = [0,0,0,0,csy+2*cym,0,-cym,-cym;
%     0,0,0,0,0,l2^2*csx+lm^2*cym,-lm*cym,lm*cym;
%     0,0,0,0,-cym,-lm*cym,cym,0;
%     0,0,0,0,-cym,lm*cym,0,cym ];
% v = 100;                        % Velocity in m/s
% w = 1/v;
% CA = [zeros(4,8);
%     C2];
% CD = [2*fy*w,0,0,0,0,0,0,0;
%     0,2*l0^2*fx*w,0,0,0,0,0,0;
%     0,0,2*fy*w,0,0,0,0,0;
%     0,0,0,2*l0^2*fx*w,0,0,0,0;
%     zeros(4,8)];
% le = 1;
% E = [0 0; 0 0; 0 0; 0 0; 1 1; le -le; 0 0; 0 0];
% A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
% 
% % With initial condition
% % Set initial condition x_o(0)
% a = 0.05;                                             % Initial condition 
% x_o = a*ones(16,1);                      
% 
% % Create new matrix B
% B = [zeros(8,2);M\E];
% Bo = B;
% % z0 = norm(x_o)/norm(B(:,1));
% z0 = a;
% B = [B,x_o/z0];                                        % x_o/z0 is a unit vector
% 
% % Controllability check
% p = 4;
% q = 12;
% C = [eye(p),zeros(p,q)];
% D = zeros(p,size(B,2));
% sys = ss(A,B,C,D);
% 
% [sysb,g,T,Ti] = balreal(sys);
% 
% order_s = 5;                                            % Truncated order
% order_e = 16;
% % sysmdc = modred(sysb,order_s:order_e,'MathDC');       % Reduced order subsystem
% sysdel = modred(sysb,order_s:order_e,'Truncate');       % Reduced order subsystem
% order = order_s-1;
% 
% Ti = Ti(1:order,:);
% T = T(:,1:order);
% 
% % Reduced system
% A_r = sysdel.A;
% B_r = sysdel.B;
% B_r = B_r(:,1:2);
% C_r = sysdel.C;
% x_r = Ti*x_o;
% 
% % Feedback gain
% Q = 10*eye(order);
% R = [10,0;0,10]/10^7;
% [K1,S1,P1] = lqr(A_r,B_r,Q,R);
% 
% ite = 0;
% for lambda = 0.1:0.01:0.2
%     ite = ite + 1;
%     kgy=lambda*W/l0;
%     Ks = [ kpy+kgy,0,0,0,-kpy,-b*kpy,0,0;
%            0,l1^2*kpx+kgfai,0,0,0,-l1^2*kpx,0,0;
%            0,0,kpy+kgy,0,-kpy,b*kpy,0,0;
%            0,0,0,l1^2*kpx+kgfai,0,-l1^2*kpx,0,0;
%            -kpy,0,-kpy,0,(2*kpy+ksy+2*kym),0,-kym,-kym;
%            -b*kpy,-l1^2*kpx,b*kpy,-l1^2*kpx,0,l2^2*ksx+2*l1^2*kpx+2*b^2*kpy+2*lm^2*kym,-lm*kym,lm*kym;
%            0,0,0,0,-kym,-lm*kym,kym,0;   
%            0,0,0,0,-kym,lm*kym,0,kym ];
%  
%     Kus = [ 0,-2*fy,0,0,0,0,0,0;
%         2*lambda*l0*fx/r0,0,0,0,0,0,0,0;
%         0,0,0,-2*fy,0,0,0,0;
%         0,0,2*lambda*l0*fx/r0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0;   
%         0,0,0,0,0,0,0,0]; 
%     A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
%     
%     eig_A(ite) = max(real(eig(A-Bo*(K1*pinv(C_r)*C))));
%     lambda_var(ite) = lambda;
% end
% 
% swFigSize
% swEPSfigure
% plot(lambda_var, eig_A)
% xlabel('$\lambda_e$')
% ylabel('$\zeta_{max}$')
% title('Stable case')
% set(gcf, 'Position', [4.874200000000000e+03 2.542000000000000e+02 560 420])