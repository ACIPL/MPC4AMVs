%% ------------------------------------------------------------------------
% License information: You are free to use or extend the code for research
% or education purposes provided that you provide clear a reference
% including:
% [1] Shen, Chao, Brad Buckham, and Yang Shi. "Modified C/GMRES algorithm for 
%     fast nonlinear model predictive tracking control of AUVs." 
%     IEEE Transactions on Control Systems Technology 25.5 (2016): 1896-1904.
% [2] Shi, Y., Shen, C., Wei, H., Zhang, K., 2022. Advanced model
%     predictive control for autonomous marine vehicels, Springer.
%
% Code description: 
% Simulation for Chapter 8 
%% ------------------------------------------------------------------------
addpath('./ch8header02/')
clc
clear all

N=5;  % prediction horizon N
T=0.1; % sampling period T
nx=6; % X=[x;y;psi;u;v;r]
nu=3; % Ui=[Fu;Fv;Fr]

Fu_max=2000;
Fv_max=2000;
Fr_max=900;
Umax=[Fu_max;Fv_max;Fr_max];

q11=1e4; q22=1e4; q33=10; q44=10; q55=1; q66=10;
Q=diag([q11,q22,q33,q44,q55,q66]);

r11=1e-4; r22=1e-4; r33=1e-2;
R=diag([r11,r22,r33]);

tau=1e2;

qf11=10; qf22=10; qf33=1; qf44=1; qf55=1; qf66=1;
Qf=diag([qf11,qf22,qf33,qf44,qf55,qf66]);
weights = {Q, R, Qf};

m=116; Iz=13.1; X_udot=-167.6; Y_vdot=-477.2; N_rdot=-15.9;
Xu=26.9; Yv=35.8; Nr=3.5; Du=241.3; Dv=503.8; Dr=76.9;

Mx=m-X_udot; My=m-Y_vdot; Mpsi=Iz-N_rdot;
disturbance=[0;0;0];

coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr];
ndof = 3;
X0=[0.5;0;0;0;0;0];
U0 = [0;0;0];
auv1 = AUV(coef, ndof, X0, U0);
auv1_kinematic = AUV_kinematic(auv1,weights);

%==========================================================================
% Tracking Reference
%==========================================================================
% ti=0:0.1:20;% size of 101
% xd=0.5*ti;% size of 101
% yd=sin(xd);  % size of 101

% ti=0:0.1:15;
% xd=0.8*cos(0.5*ti);
% yd=0.8*sin(0.5*ti);

ti=0:0.1:30;
xd=sin(ti/2);
yd=sin(ti/4);

psi_d=atan2(diff(yd),diff(xd));
xd=xd';
yd=yd';
psi_d=psi_d';

xd_dot=diff(xd)/T;
yd_dot=diff(yd)/T;
xd_2dot=diff(xd_dot)/T;
yd_2dot=diff(yd_dot)/T;
L_2dot=length(xd_2dot);
ud=zeros(L_2dot,1);

for i=1:1:L_2dot
    ud(i)=sqrt(xd_dot(i)^2+yd_dot(i)^2);
end

vd=zeros(L_2dot,1);
rd=zeros(L_2dot,1);

for i=1:1:L_2dot
    rd(i)=(xd_dot(i)*yd_2dot(i)-yd_dot(i)*xd_2dot(i))/(xd_dot(i)^2+yd_dot(i)^2);
end

M=L_2dot-N;
Pall=[xd(1:L_2dot),yd(1:L_2dot),psi_d(1:L_2dot),ud,vd,rd]';

%==========================================================================
X=zeros(nx,N+1); % X=[X0,X1,...,XN]
U=zeros(nu,N);   % U=[U0,U1,...,U_N-1]

Xall=zeros(nx,M+1);
Uall=zeros(nu,M);
X0=[0.3;0;0;0;0;0];
Xall(:,1)=X0;
U0=Uall(:,1:N);

iter_barrier=6;
for i=1:1:iter_barrier
    P=Pall(:,1:1+N);
    tic
    U1 = fsolve(@(U) auv1_kinematic.F_AUV_kinematic( X0,U,P,N,T,Umax,tau),U0);
    toc
    U0=U1;
    tau=tau/10;
    i
end

tic
for i=1:1:M
    tau=1e2;
    P=Pall(:,i:i+N);
    
    for j=1:1:3
        U1 = fsolve(@(U) auv1_kinematic.F_AUV_kinematic( X0,U,P,N,T,Umax,tau),U0);
        U0=U1;
        tau=tau/10;
    end
    Uall(:,i)=U1(:,1);
    auv1.advance(Uall(:,i),disturbance,T);
    Xplus = auv1.X;
    Xall(:,i+1)=Xplus;
    U0=[U1(:,2:N),zeros(nu,1)];
    X0=Xplus;
    i
end
toc

xd2=xd';
yd2=yd';

xe=Xall(1,1:M)-xd2(1,1:M);
ye=Xall(2,1:M)-yd2(1,1:M);

%==========================================================================
% plot
figure(1)
plot(Xall(1,:),Xall(2,:),'k');
grid on;

figure(2)
subplot(3,1,1)
plot(ti(1,1:M),Uall(1,:),'k'), title('Fu')
grid on;
subplot(3,1,2)
plot(ti(1,1:M),Uall(2,:),'k'), title('Fv')
grid on;
subplot(3,1,3)
plot(ti(1,1:M),Uall(3,:),'k'), title('Fr')
grid on;

figure(3)
subplot(2,1,1)
plot(ti(1,1:M),xe(1,1:M),'k'), title('xe')
grid on;
subplot(2,1,2)
plot(ti(1,1:M),ye(1,1:M),'k'), title('ye')
grid on;

rmpath('./ch8header02/')

