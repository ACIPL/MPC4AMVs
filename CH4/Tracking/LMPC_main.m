%% ------------------------------------------------------------------------
% License information: You are free to use or extend the code for research
% or education purposes provided that you provide clear a reference
% including:
% [1] Shen, Chao, Yang Shi, and Brad Buckham. "Trajectory tracking control 
%     of an autonomous underwater vehicle using Lyapunov-based model 
%     predictive control." IEEE Transactions on Industrial Electronics 65.7 (2017): 5796-5805.
%
% Code description: Simulation of the MPC-based Trajectory Tracking.
% -------------------------------------------------------------------------
addpath('./utilis/','./ch4header/')
clc;
clear;

%======================================================================%
% System coefficients
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

B=[0.7974    0.8643    0.8127    0.8270
    0.6032    0.5029   -0.5824   -0.5610
    0.2945   -0.3302   -0.2847    0.3505];
Bplus=pinv(B);

M=diag([Mx;My;Mpsi]);
PI=[0,-1,0;1,0,0;0,0,0];
g=0;
Lambda=diag([1;1;1]);

% disturbance=[10;10;0];
disturbance=[0;0;0];

%======================================================================%
% create AUV objects
%======================================================================%
coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr];
ndof = 3;
X0=[0.5;0;0;0;0;0];
Thrust0 = [0;0;0;0];
auv1 = AUV(coef, ndof, B, X0, Thrust0);
auv2 = AUV(coef, ndof, B, X0, Thrust0);

%======================================================================%
% create a nonlinear BS controller object
%======================================================================%
Kp=diag([1;1;1]);
Kd=diag([1;1;1]);
bs1 = nonlinearBS_controller(Kp,Kd);

%======================================================================%
% create a sine2D object
%======================================================================%
% Case 1: a sinusoidal reference trajectory
ref_trajectory = sine2D(1,0.5,0); 
% Case 2: an eight-shaped reference trajectory
% ref_trajectory = eight_shape2D(1,0.5,0.25,0); 
%======================================================================%
% Nonlinear BS control simulation
%======================================================================%

ref = [ref_trajectory.xR;ref_trajectory.yR;ref_trajectory.psiR;...
    ref_trajectory.xRd;ref_trajectory.yRd;ref_trajectory.psiRd;...
    ref_trajectory.xRdd;ref_trajectory.yRdd;ref_trajectory.psiRdd];


dt=0.1;
Tstep=200;
nx=length(X0);
nu=3;
n_thrust=4;

Xall=zeros(nx,Tstep+1);
Uall=zeros(nu,Tstep);
Thrust_all=zeros(n_thrust,Tstep);

eta_Ref_all=zeros(9,Tstep);

Xall(:,1)=X0;
Xplus=X0;

for i=1:1:Tstep
    eta_Ref_all(:,i) = ref;
    
    tau = bs1.calc_control(ref,auv1.X);
    Uall(:,i)= tau;
    thrust = Bplus*tau;
    Thrust_all(:,i) = thrust;
    auv1.advance(thrust,disturbance,dt);

    Xall(:,i+1) = auv1.X;
    ref_trajectory.update(dt);
    ref_trajectory.t
    ref = [ref_trajectory.xR;ref_trajectory.yR;ref_trajectory.psiR;...
        ref_trajectory.xRd;ref_trajectory.yRd;ref_trajectory.psiRd;...
        ref_trajectory.xRdd;ref_trajectory.yRdd;ref_trajectory.psiRdd];
    
end
figure(1)
plot(eta_Ref_all(1,:),eta_Ref_all(2,:),'k');
hold on;
plot(Xall(1,:),Xall(2,:),'b');
hold on;
%======================================================================%
% create a LMPC tracking controller object
%======================================================================%
% Some intial setup for LMPC

N=5;  % prediction horizon N
Tstep2 = Tstep-N;

% Weighting Matrices
q11=1e5;
q22=1e5;
q33=1e3;
q44=1e2;
q55=1e2;
q66=1e2;
Q=diag([q11,q22,q33,q44,q55,q66]);

r11=1e-4;
r22=1e-4;
r33=1e-4;
r44=1e-4;
R=diag([r11,r22,r33,r44]);

R1=0;

qf11=1e3;
qf22=1e3;
qf33=1e2;
qf44=1e1;
qf55=1e1;
qf66=1e1;
Qf=diag([qf11,qf22,qf33,qf44,qf55,qf66]);

weights = {Q,R,Qf};

T1_max = 500;
T2_max = 500;
T3_max = 500;
T4_max = 500;
Thrust_max0 = [T1_max;T2_max;T3_max;T4_max];
Thrust_min0 = [-T1_max;-T2_max;-T3_max;-T4_max];

internal_model = AUV(coef, ndof, B, X0, Thrust0);
auxiliary_controller = nonlinearBS_controller(Kp,Kd);
lmpc1 = LMPC_controller(N,internal_model,auxiliary_controller,weights,Thrust_max0,Thrust_min0);

%======================================================================%
% LMPC simulation
%======================================================================%

Xall2 = zeros(nx,Tstep2+1);
Uall2 = zeros(nu,Tstep2);
Thrust_all2 = zeros(n_thrust,Tstep);
Xall2(:,1) = X0;
U0 = Uall2(:,1);
Thrust0 = Bplus*U0;

M = nx + nu + n_thrust;
t = 0;

u0 = zeros(n_thrust*N,1);

tic
for i=1:1:Tstep2
        
    t1 = t;
    TIme(i,1) = t;
    
    P = RefGen(M,N,t1,dt,coef,Bplus,ref_trajectory);
    u = lmpc1.calc_control(P,X0,B,u0,dt,ref_trajectory,t1);
    
    u_actual = u(1:n_thrust,1);
    Thrust_all2(:,i) = u_actual;
    
    tau_actual = B*u_actual;
    Uall2(:,i) = tau_actual;
   
    auv2.advance(u_actual,disturbance,dt);
    Xall2(:,i+1) = auv2.X;
    X0 = auv2.X;
    %======================================================================%
    % Construction of a feasible control u0_fmincon using PD controller 
    %======================================================================%
    t2 = t+dt;
    P_ = RefGen(M,N,t2,dt,coef,Bplus,ref_trajectory);
    u0 = lmpc1.calc_initial_guess(P_,X0,Bplus,dt);  
    i
    t = t + dt;
end

%==========================================================================
% Plot

figure(1)
plot(eta_Ref_all(1,:),eta_Ref_all(2,:),'k');
hold on;
plot(Xall(1,1:Tstep2),Xall(2,1:Tstep2),'b');
hold on;
plot(Xall2(1,1:Tstep2),Xall2(2,1:Tstep2),'r');
grid on;
legend('Ref', 'Nonlinear BS Control', 'LMPC-Based Control');
xlabel('x [m]');
ylabel('y [m]');

figure(2)
subplot(3,1,1)
plot(TIme(1:Tstep2,1),Uall(1,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Uall2(1,1:Tstep2),'r');
grid on;
ylabel('F_u [N]');
legend('BS','LMPC');
subplot(3,1,2)
plot(TIme(1:Tstep2,1),Uall(2,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Uall2(2,1:Tstep2),'r');
grid on;
ylabel('F_v [N]');

subplot(3,1,3)
plot(TIme(1:Tstep2,1),Uall(3,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Uall2(3,1:Tstep2),'r');
grid on;
xlabel('Time [s]');
ylabel('F_r [Nm]');


figure(3)
subplot(4,1,1)
plot(TIme(1:Tstep2,1),Thrust_all(1,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(1,1:Tstep2),'r');
grid on;
ylabel('u_1 [N]')

subplot(4,1,2)
plot(TIme(1:Tstep2,1),Thrust_all(2,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(2,1:Tstep2),'r');
grid on;
ylabel('u_2 [N]')

subplot(4,1,3)
plot(TIme(1:Tstep2,1),Thrust_all(3,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(3,1:Tstep2),'r');
grid on;
ylabel('u_3 [N]')

subplot(4,1,4)
plot(TIme(1:Tstep2,1),Thrust_all(4,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(4,1:Tstep2),'r');
grid on;
ylabel('u_4 [N]')
xlabel('Time [s]');

figure(4)
subplot(3,1,1)
plot(TIme(1:Tstep2,1),eta_Ref_all(1,1:Tstep2),'k');
hold on;
plot(TIme(1:Tstep2,1),Xall(1,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(1,1:Tstep2),'r');
ylabel('x [m]');

subplot(3,1,2)
plot(TIme(1:Tstep2,1),eta_Ref_all(2,1:Tstep2),'k');
hold on;
plot(TIme(1:Tstep2,1),Xall(2,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(2,1:Tstep2),'r');
ylabel('y [m]');

subplot(3,1,3)
plot(TIme(1:Tstep2,1),eta_Ref_all(3,1:Tstep2),'k');
hold on;
plot(TIme(1:Tstep2,1),Xall(3,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(3,1:Tstep2),'r');
xlabel('Time [s]');
ylabel('\psi [rad]');


% figure(5)
% subplot(3,1,1)
% plot(TIme(1:Tstep2,1),Xall(4,1:Tstep2),'b');
% hold on;
% grid on;
% plot(TIme(1:Tstep2,1),Xall2(4,1:Tstep2),'r');
% 
% 
% subplot(3,1,2)
% plot(TIme(1:Tstep2,1),Xall(5,1:Tstep2),'b');
% hold on;
% grid on;
% plot(TIme(1:Tstep2,1),Xall2(5,1:Tstep2),'r');
% 
% 
% subplot(3,1,3)
% plot(TIme(1:Tstep2,1),Xall(6,1:Tstep2),'b');
% hold on;
% grid on;
% plot(TIme(1:Tstep2,1),Xall2(6,1:Tstep2),'r');


MSEx1=0;
MSEy1=0;
MSEpsi1=0;

MSEx2=0;
MSEy2=0;
MSEpsi2=0;

for i=1:1:Tstep2
    
    MSEx1=MSEx1+(Xall(1,i)-eta_Ref_all(1,i))^2;
    MSEx2=MSEx2+(Xall2(1,i)-eta_Ref_all(1,i))^2;
    
    MSEy1=MSEy1+(Xall(2,i)-eta_Ref_all(2,i))^2;
    MSEy2=MSEy2+(Xall2(2,i)-eta_Ref_all(2,i))^2;
    
    MSEpsi1=MSEpsi1+(Xall(3,i)-eta_Ref_all(3,i))^2;
    MSEpsi2=MSEpsi2+(Xall2(3,i)-eta_Ref_all(3,i))^2;
    
end

MSEx1=MSEx1/Tstep2;
MSEx2=MSEx2/Tstep2;

MSEy1=MSEy1/Tstep2;
MSEy2=MSEy2/Tstep2;

MSEpsi1=MSEpsi1/Tstep2;
MSEpsi2=MSEpsi2/Tstep2;

rmpath('./utilis/','./ch4header/')