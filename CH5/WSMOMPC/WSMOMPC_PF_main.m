%% ------------------------------------------------------------------------
% License information: You are free to use or extend the code for research
% or education purposes provided that you provide clear a reference
% including:
% [1] Shen, Chao, Yang Shi, and Brad Buckham. "Path-following control of an 
%     AUV: A multiobjective model predictive control approach." IEEE
%     Transactions on Control Systems Technology 27.3 (2018): 1334-1342.
% [2] Shi, Y., Shen, C., Wei, H., Zhang, K., 2022. Advanced model
%     predictive control for autonomous marine vehicels, Springer.
%
% Code description: 
% Simulation for Section 5.5.2 WS-MOMPC Design for trajectory tracking control
%% ------------------------------------------------------------------------
addpath('./utilis/','./ch5WSheader/')
clc;
clear;

%======================================================================%
% System coefficients & Other parameters
%======================================================================%
err_model=0;
m=116;
Iz=13.1;
X_udot=-167.6*(1+err_model);
Y_vdot=-477.2*(1+err_model);
N_rdot=-15.9*(1+err_model);
Xu=26.9*(1+err_model);
Yv=35.8*(1+err_model);
Nr=3.5*(1+err_model);
Du=241.3*(1+err_model);
Dv=503.8*(1+err_model);
Dr=76.9*(1+err_model);

Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;

% disturbance=[10;10;0];
disturbance=[0;0;0;0];

X0=[0.5;0;0;0;0;0];
s0=0;
dt=0.1;
Tstep=100;
nx=length(X0);
nu=4; % Ui=[Fu;Fv;Fr,vs]
uR=1;
N=10;  % prediction horizon N

Xall=zeros(nx,Tstep+1);
Uall=zeros(nu,Tstep);
Sall=zeros(1,Tstep+1);
U0=Uall(:,1);

Xall(:,1)=X0;
Xplus=X0;

alpha=0.99;
beta=1e0;

s_plus=s0;
xR=s0;
yR=sin(s0);
psiR=atan(cos(s0));
%======================================================================%
% create AUV objects
%======================================================================%
coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr];
ndof = 3;
auv1 = AUVI(coef, ndof, X0, U0);

% Weighting Matrices
Q1=diag([1e5,1e5,1e2,1e-1,1e-1,1e-1,1e-1]);
Q2 = diag([1e0,1e0,1e0,1e3,1e-1,1e-1,1e-1]);

R1=diag([1e-3,1e-3,1e-3,1e-3]);
R2=diag([1e-3,1e-3,1e-3,1e-3]);

Qf1=diag([1e2,1e2,1e1,1e-3,1e-3,1e-3,1e-3]);
Qf2 = diag([1e-1,1e-1,1e-1,1e2,1e-3,1e-3,1e-3]);
weights = {Q1,R1,Qf1,Q2,R2,Qf2};

Fu_max=500;
Fv_max=500;
Fr_max=500;
v_smax=100;
v_smin=-100;
Umax0=[Fu_max;Fv_max;Fr_max;v_smax];
Umin0=[-Fu_max;-Fv_max;-Fr_max;v_smin];

internal_model = AUVI(coef, ndof, X0, U0);
path = sine2D(1,1,0);
wsmompc = MPC_controller(N,internal_model,weights,Umax0,Umin0,path);

%======================================================================%
% WSMOMPC simulation
%======================================================================
u0 = zeros(nu*N,1);
tic
for i=1:1:Tstep
    err=(Xall(1,i)-xR)^2+(Xall(2,i)-yR)^2+(Xall(3,i)-psiR)^2;
    alpha=alpha_logistic(err, beta);
    
    u = wsmompc.calc_control(uR,alpha,s_plus,X0,u0,dt);
    u_actual = u(1:nu,1);
    Uall(:,i) = u_actual;
   
    auv1.advance(u_actual,disturbance,dt);
    Xall(:,i+1) = auv1.X;
    X0 = auv1.X;
    s_plus = path_evol(s_plus,u(4,1),dt);
    Sall(i+1)=s_plus;
    %=====================================================================%
    % Construction of a feasible control u0
    %=====================================================================%
    u0 = [u(nu+1:nu*N,1);u(nu*(N-1)+1:nu*N,1)];  
    xR=s_plus;
    yR=sin(s_plus);
    psi_R=atan(cos(s_plus));
    i
end

%==========================================================================
% Plot
xs=Sall;
ys=Sall;
for i=1:1:Tstep+1
    
    ys(i)=sin(xs(i));
    
end

OT=zeros(Tstep,1);

for i=1:1:Tstep
    OT(i)=i*dt;
end

% plot figures
figure(1)
xd=0:0.1:10; % size of 101
yd=sin(xd);  % size of 101
Pall=[xd;yd]';
plot(Pall(:,1),Pall(:,2),'g')
hold on;
plot(xs,ys,'r')
hold on;
plot(Xall(1,:),Xall(2,:),'b');
hold on;

figure(2)
subplot(2,2,1)
plot(OT(:,1),Uall(1,:),'b'), title('Fu')
hold on;
grid on;

subplot(2,2,2)
plot(OT(:,1),Uall(2,:),'b'), title('Fv')
hold on;
grid on;

subplot(2,2,3)
plot(OT(:,1),Uall(3,:),'b'), title('Fr')
hold on;
grid on;

subplot(2,2,4)
plot(OT(:,1),Uall(4,:),'b'), title('vs')
hold on;
grid on;

figure(3)
plot(Xall(1,:),Xall(2,:),'b');
hold on;
grid on;

figure(4)
plot(OT(:,1),Xall(4,1:Tstep));
title('Surge Velocity of the Vehicle')
ylabel('u')
hold on;
grid on;

rmpath('./utilis/','./ch5WSheader/')