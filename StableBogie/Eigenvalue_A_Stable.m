%=============================================================%
% Eigenvalue changing of the system
% Robustness check
%=============================================================%

clc
clear
close all
format long
% swEPSfigure

W= 21000*9.81;
mm=3530/2;     % 1765kg                  
Im=1019/2;

mw=2585;
Iw=2024;
lm=0.5;

mf=3656+mm*2;               
If=7849+2*mm*lm^2+Im*2;        
    
lambda=0.1;
csx=1000000*1;
csy=60000*1;
kpx=4*10^7;  % 4*10^7;
l0=0.7465;
kpy=6*10^6; %+lambda*W/l0;
kgy=lambda*W/l0;

kgfai=-23000;       
ksx=0.2*10^6*4 ;   
ksy=0.2*10^6*4 ;    
fy=8.62400*10^6;
fx=8.14400*10^6;

l1=1.1;
l2=0.95;
b=2.5/2;
r0=0.625;

freqmy=1.8;
dumpmy=0.2;

kym=(2*pi*freqmy)^2*mm;
cym=dumpmy*2*mm*(2*pi*freqmy);
M0=[mw,Iw,mw,Iw,mf,If,mm,mm];
M=diag(M0);

Ks = [ kpy+kgy,0,0,0,-kpy,-b*kpy,0,0;
     0,l1^2*kpx+kgfai,0,0,0,-l1^2*kpx,0,0;
     0,0,kpy+kgy,0,-kpy,b*kpy,0,0;
     0,0,0,l1^2*kpx+kgfai,0,-l1^2*kpx,0,0;
     -kpy,0,-kpy,0,(2*kpy+ksy+2*kym),0,-kym,-kym;
     -b*kpy,-l1^2*kpx,b*kpy,-l1^2*kpx,0,l2^2*ksx+2*l1^2*kpx+2*b^2*kpy+2*lm^2*kym,-lm*kym,lm*kym;
     0,0,0,0,-kym,-lm*kym,kym,0;   
     0,0,0,0,-kym,lm*kym,0,kym ];
 
Kus = [ 0,-2*fy,0,0,0,0,0,0;
     2*lambda*l0*fx/r0,0,0,0,0,0,0,0;
     0,0,0,-2*fy,0,0,0,0;
     0,0,2*lambda*l0*fx/r0,0,0,0,0,0;
     0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0;   
     0,0,0,0,0,0,0,0]; 
 
C2 = [0,0,0,0,csy+2*cym,0,-cym,-cym;
    0,0,0,0,0,l2^2*csx+lm^2*cym,-lm*cym,lm*cym;
    0,0,0,0,-cym,-lm*cym,cym,0;
    0,0,0,0,-cym,lm*cym,0,cym ];
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

% With initial condition
% Set initial condition x_o(0)
a = 0.05;                                             % Initial condition 
x_o = a*ones(16,1);                      

% Create new matrix B
B = [zeros(8,2);M\E];
Bo = B;
% z0 = norm(x_o)/norm(B(:,1));
z0 = a;
B = [B,x_o/z0];                                        % x_o/z0 is a unit vector

% Balanced system
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

Ti = Ti(1:order,:);
T = T(:,1:order);

% Reduced system
A_r = sysdel.A;
B_r = sysdel.B;
B_r = B_r(:,1:2);
C_r = sysdel.C;
x_r = Ti*x_o;

% Feedback gain
Q = 100000*eye(order);
R = [10,0;0,10]/10000000;
[K1,S1,P1] = lqr(A_r,B_r,Q,R);

ite = 0;
for v = 50:150
    ite = ite + 1;
    w = 1/v;
    CD = [2*fy*w,0,0,0,0,0,0,0;
        0,2*l0^2*fx*w,0,0,0,0,0,0;
        0,0,2*fy*w,0,0,0,0,0;
        0,0,0,2*l0^2*fx*w,0,0,0,0;
        zeros(4,8)];
    A = [zeros(8),eye(8);-inv(M)*(Ks+Kus),-inv(M)*(CA+CD)];
    
    eig_A(ite) = max(real(eig(A-Bo*(K1*pinv(C_r)*C))));
    v_var(ite) = v*3.6;
end

swFigSize
swEPSfigure
plot(v_var, eig_A)
xlabel('Velcity (km/h)')
ylabel('Max(real(eig(A)))')

% print -depsc EigAChanging_Stablecase.eps
