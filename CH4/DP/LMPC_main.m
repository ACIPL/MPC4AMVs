%% ------------------------------------------------------------------------
% License information: You are free to use or extend the code for research
% or education purposes provided that you provide clear a reference
% including:
% [1] Shen, Chao, Yang Shi, and Brad Buckham. "Lyapunov-based model predictive 
% control for dynamic positioning of autonomous underwater vehicles." 
% 2017 IEEE International Conference on Unmanned Systems (ICUS). IEEE, 2017.
%
% Code description: Simulation of the MPC-based DP.
% -------------------------------------------------------------------------
addpath('./utils/','./ch4DPheader/')

clc;
clear;
%======================================================================%
% System coefficients
%======================================================================%
err_model = 0; % err_model=0.2 implies AUV is subject to 20% model error
m = 116;
Iz = 13.1;
X_udot = -167.6 * (1 + err_model);
Y_vdot = -477.2 * (1 + err_model);
N_rdot = -15.9 * (1 + err_model);
Xu = 26.9 * (1 + err_model);
Yv = 35.8 * (1 + err_model);
Nr = 3.5 * (1 + err_model);
Du = 241.3 * (1 + err_model);
Dv = 503.8 * (1 + err_model);
Dr = 76.9* (1 + err_model);
Mudot = m - X_udot;
Mvdot = m - Y_vdot;
Mrdot = Iz - N_rdot;

B = [0.7974    0.8643    0.8127    0.8270
0.6032    0.5029   -0.5824   -0.5610
0.2945   -0.3302   -0.2847    0.3505];

Bplus = pinv(B);
dt = 0.1;
disturbance = [0;0;0]; % disturbances experting on AUV

%======================================================================%
% create AUV objects
%======================================================================%
coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr];
ndof = 3;
X0 = [5;5;-pi/2;0;0;0];
Thrust0 = [0;0;0;0];
auv1 = AUV(coef, ndof, B, X0, Thrust0);
auv2 = AUV(coef, ndof, B, X0, Thrust0);
%======================================================================%
% create a nonlinear PD controller object
%======================================================================%
Kp = diag([10;10;10]);
Kd = diag([10;10;10]);
pd1 = nonlinearPD_controller(Kp,Kd);

%======================================================================%
% create a point2D object
%======================================================================%
ref_trajectory = point2D(0,0,0); % (xd,yd,psid)

%======================================================================%
% Nonlinear PD control simulation
%======================================================================%
ref = [ref_trajectory.xR;ref_trajectory.yR;ref_trajectory.psiR];

Tstep = 700;
nx = length(X0);
nu = 3;
n_thrust = 4;
Xall = zeros(nx,Tstep+1);
Uall = zeros(nu,Tstep);
Thrust_all = zeros(n_thrust,Tstep);
Xall(:,1) = X0;

for i=1:1:Tstep
    tau = pd1.calc_control(ref,auv1.X);  
    Uall(:,i)= tau;
    thrust = Bplus*tau;
    Thrust_all(:,i) = thrust;   
    auv1.advance(thrust,disturbance,dt);
    Xall(:,i+1) = auv1.X;    
    ref_trajectory.update(dt);
    ref = [ref_trajectory.xR;ref_trajectory.yR;ref_trajectory.psiR];
end

%======================================================================%
% create a LMPC controller object
%======================================================================%
N = 3;
Tstep2 = Tstep - N; 
TIme = zeros(Tstep2,1);

q11=1e5;
q22=1e5;
q33=1e4;
q44=1e3;
q55=1e3;
q66=1e3;
Q=diag([q11,q22,q33,q44,q55,q66]);
r11=1e-3;
r22=1e-3;
r33=1e-3;
r44=1e-3;
R=diag([r11,r22,r33,r44]);
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
auxiliary_controller = nonlinearPD_controller(Kp,Kd);
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
    u = lmpc1.calc_control(P,X0,B,u0,dt);
    
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
toc
%======================================================================%
% Plot figures
%======================================================================%
% The AUV trajectories
figure(1)
plot(Xall(1,1:Tstep),Xall(2,1:Tstep),'b');
hold on;
plot(Xall2(1,1:Tstep2),Xall2(2,1:Tstep2),'r');
grid on;
legend('Nonlinear PD Control','LMPC-Based Control','Location','NorthWest');
xlabel('x [m]');
ylabel('y [m]');
auv1.X
auv2.X

% The state trajectories
figure(2)
subplot(3,1,1)
plot(0:dt:dt*Tstep,Xall(1,:),'b');
hold on;
plot(0:dt:dt*Tstep2,Xall2(1,:),'r');
grid on;
legend('Nonlinear PD Control','LMPC-Based Control','Location','NorthEast');
ylabel('x [m]');
subplot(3,1,2)
plot(0:dt:dt*Tstep,Xall(2,:),'b');
hold on;
plot(0:dt:dt*Tstep2,Xall2(2,:),'r');
ylabel('y [m]');
grid on;
subplot(3,1,3)
plot(0:dt:dt*Tstep,Xall(3,:),'b');
hold on;
grid on;
plot(0:dt:dt*Tstep2,Xall2(3,:),'r');
xlabel('Time [s]');
ylabel('\psi [rad]');

% The state trajectories
figure(3)
subplot(4,1,1)
plot(0:dt:dt*(Tstep-1),Thrust_all(1,:),'b');
hold on;
plot(0:dt:dt*(Tstep-1),Thrust_all2(1,:),'r');
grid on;
legend('Nonlinear PD Control','LMPC-Based Control','Location','NorthEast');
ylabel('u_1 [N]');
subplot(4,1,2)
plot(0:dt:dt*(Tstep-1),Thrust_all(2,:),'b');
hold on;
plot(0:dt:dt*(Tstep-1),Thrust_all2(2,:),'r');
grid on;
ylabel('u_2 [N]');
subplot(4,1,3)
plot(0:dt:dt*(Tstep-1),Thrust_all(3,:),'b');
hold on;
plot(0:dt:dt*(Tstep-1),Thrust_all2(3,:),'r');
grid on;
ylabel('u_3 [N]');
subplot(4,1,4)
plot(0:dt:dt*(Tstep-1),Thrust_all(3,:),'b');
hold on;
plot(0:dt:dt*(Tstep-1),Thrust_all2(3,:),'r');
grid on;
xlabel('Time [s]');
ylabel('u_4 [N]');

% delete the added path
rmpath('./utils/','./ch4DPheader/')
















