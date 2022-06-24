% Simulation for Section 6.4.2 DLMPC for formation tracking control
% First scenario: Formation tracking without considering the collision
% avoidance, i.e., the collision avoidance cost term is not involved in the
% overall objective function
% 2021/06/29 @ ELW 133, UVic

addpath('./utilis/','./ch6header/')
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

M=diag([Mx;My;Mpsi]);
PI=[0,-1,0;1,0,0;0,0,0];
g=0;
Lambda=diag([1;1;1]);
disturbance=[0;0;0];

%======================================================================%
% create 3 AUVs
%======================================================================%
coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr];
ndof = 3;
X1_0 = [0; 2; pi/2; 0.4; 0; 0];
X2_0 = [-0.5; -0.5; 3*pi/4; 0.6; 0; 0];
X3_0 = [-0.1; 0.1; -pi/2; 0.7; 0; 0];
U0 = [0;0;0];
% Denote the actual system
auv1 = AUV(coef, ndof, X1_0, U0);
auv2 = AUV(coef, ndof, X2_0, U0);
auv3 = AUV(coef, ndof, X3_0, U0);
% Denote the prediction model
internal_auv1 = AUV(coef, ndof, X1_0, U0);
internal_auv2 = AUV(coef, ndof, X2_0, U0);
internal_auv3 = AUV(coef, ndof, X3_0, U0);

%======================================================================%
% create 3 nonlinear BS controller instances
%======================================================================%
Kp=diag([1;1;1]);
Kd=diag([1;1;1]);
bs1 = nonlinearBS_controller(Kp,Kd);
bs2 = nonlinearBS_controller(Kp,Kd);
bs3 = nonlinearBS_controller(Kp,Kd);

%======================================================================%
% create three sine2D reference trajectories
%======================================================================%
t0 = 0;
ref_trajectory = sine2D(1,0.5,pi/2, t0);
ref = [ref_trajectory.xR;ref_trajectory.yR;ref_trajectory.psiR;...
    ref_trajectory.xRd;ref_trajectory.yRd;ref_trajectory.psiRd;...
    ref_trajectory.xRdd;ref_trajectory.yRdd;ref_trajectory.psiRdd];
d1r = [0; 1; 0; 0; 0; 0; 0; 0; 0];
d2r = [0; -1; 0; 0; 0; 0; 0; 0; 0];
d3r = [-1; 0; 0; 0; 0; 0; 0; 0; 0];

ref1 = ref + d1r;
ref2 = ref + d2r;
ref3 = ref + d3r;

%==========================================================================
% Other design parameters
%==========================================================================
% sampling time
dt = 0.1;
% simulation step
Tstep = 150;
nx = length(X1_0);
nu = length(U0);

[Xa1, Xa2, Xa3] = deal(zeros(nx,Tstep+1));
[Ua1, Ua2, Ua3] = deal(zeros(nu,Tstep));

[Ref1, Ref2, Ref3, Ref] = deal(zeros(9,Tstep));

Xa1(:,1) = X1_0; Xa2(:,1) = X2_0; Xa3(:,1) = X3_0;
X1plus = X1_0; X2plus = X2_0; X3plus = X3_0;
%======================================================================%
% Nonlinear BS control simulation
%======================================================================%
for i=1:1:Tstep
    Ref(:,i) = ref;
    Ref1(:,i) = ref1;
    tau1 = bs1.calc_control(ref1, X1plus);
    Ua1(:,i) = tau1;
    X1plus = auv1.dynamics_discrete(X1plus, tau1, dt);
    Xa1(:,i+1) = X1plus;
    
    Ref2(:,i) = ref2;
    tau2 = bs2.calc_control(ref2, X2plus);
    Ua2(:,i) = tau2;
    X2plus = auv2.dynamics_discrete(X2plus, tau2, dt);
    Xa2(:,i+1) = X2plus;
    
    Ref3(:,i) = ref3;
    tau3 = bs3.calc_control(ref3, X3plus);
    Ua3(:,i) = tau3;
    X3plus = auv3.dynamics_discrete(X3plus, tau3, dt);
    Xa3(:,i+1) = X3plus;
    
    ref_trajectory.update(dt);
    ref_trajectory.t
    ref = [ref_trajectory.xR;ref_trajectory.yR;ref_trajectory.psiR;...
        ref_trajectory.xRd;ref_trajectory.yRd;ref_trajectory.psiRd;...
        ref_trajectory.xRdd;ref_trajectory.yRdd;ref_trajectory.psiRdd];
    ref1 = ref + d1r;
    ref2 = ref + d2r;
    ref3 = ref + d3r;
end


%======================================================================%
% create DLMPC tracking controllers
%======================================================================%
% Some intial setup for DLMPC
N = 4;  % prediction horizon N
M = nx + nu;
t = 0;
[u1_0, u2_0, u3_0] = deal(zeros(nu*N,1));
% Weighting Matrices
Qi = diag([1e4 1e4 1e3 1e2 1e2 1e2]);
Ri = diag([1e-3 1e-3 1e-2]);
Qfi = diag([1e4 1e4 1e3 1e2 1e2 1e2]);
weights = {Qi,Ri,Qfi};

U1_max = 1000;
U2_max = 1000;
U3_max = 1000;
U_max0 = [U1_max;U2_max;U3_max];
U_min0 = [-U1_max;-U2_max;-U3_max];

auxiliary_controller1 = nonlinearBS_controller(Kp, Kd);
auxiliary_controller2 = nonlinearBS_controller(Kp, Kd);
auxiliary_controller3 = nonlinearBS_controller(Kp, Kd);
dlmpc1 = DLMPC_controller(N,internal_auv1,auxiliary_controller1,weights,U_max0,U_min0);
dlmpc2 = DLMPC_controller(N,internal_auv2,auxiliary_controller2,weights,U_max0,U_min0);
dlmpc3 = DLMPC_controller(N,internal_auv3,auxiliary_controller3,weights,U_max0,U_min0);

%======================================================================%
% DLMPC simulation
%======================================================================%
% To store the simulation data
[X1, X2, X3] = deal(zeros(nx,Tstep+1));
[U1, U2, U3] = deal(zeros(nu,Tstep));
X1(:,1) = X1_0; X2(:,1) = X2_0; X3(:,1) = X3_0;
[X1_b, X2_b, X3_b] = deal(zeros(nx,N));
TIme = zeros(Tstep);
tic
for i=1:1:Tstep
    tic
    t1 = t;
    TIme(i+1) = t1;
    % AUV 1
    P = RefGen(M, N, t1, dt, coef, ref_trajectory);
    P1 = P + d1r;
    X1_neighbors = {X2_b, X3_b};
    [u1, X1_b] = dlmpc1.calc_control(P1, ref_trajectory, X1_neighbors, X1_0, u1_0, dt, t1, d1r);
    
    u1_actual = u1(1:nu, 1);
    U1(:,i) = u1_actual;
    
    auv1.advance(u1_actual, disturbance, dt);
    X1(:,i+1) = auv1.X;
    X1_0 = auv1.X;
    % Construction of a feasible control u1_0
    t2 = t + dt;
    P_ = RefGen(M, N, t2, dt, coef, ref_trajectory);
    P1_ = P_ + d1r;
    u1_0 = dlmpc1.calc_initial_guess(P1_, X1_0, dt);
    
    %======================================================================
    % AUV 2
    %======================================================================
    P2 = P + d2r;
    X2_neighbors = {X1_b, X3_b};
    [u2, X2_b] = dlmpc2.calc_control(P2, ref_trajectory, X2_neighbors, X2_0, u2_0, dt, t1, d2r);
    
    u2_actual = u2(1:nu, 1);
    U2(:,i) = u2_actual;
    
    auv2.advance(u2_actual, disturbance, dt);
    X2(:,i+1) = auv2.X;
    X2_0 = auv2.X;
    %======================================================================%
    % Construction of a feasible control u2_0 using BS controller
    P2_ = P_ + d2r;
    u2_0 = dlmpc1.calc_initial_guess(P2_, X2_0, dt);
    
    %======================================================================
    % AUV 3
    %======================================================================
    P3 = P + d3r;
    X3_neighbors = {X1_b, X2_b};
    [u3, X3_b] = dlmpc3.calc_control(P3, ref_trajectory, X3_neighbors, X3_0, u3_0, dt, t1, d3r);
    
    u3_actual = u3(1:nu, 1);
    U3(:,i) = u3_actual;
    
    auv3.advance(u3_actual, disturbance, dt);
    X3(:,i+1) = auv3.X;
    X3_0 = auv3.X;
    %======================================================================%
    % Construction of a feasible control u1_0 using BS controller
    P3_ = P_ + d3r;
    u3_0 = dlmpc3.calc_initial_guess(P3_, X3_0, dt);
    i
    t = t + dt;
    toc
end

%==========================================================================
% Plot

figure(1);
a1=plot(Ref(1,:),Ref(2,:), 'k-','LineWidth',1);
hold on;
a2=plot(X1(1,:), X1(2,:), '-','LineWidth', 2);
hold on;
a3=plot(X2(1,:), X2(2,:), '-','LineWidth', 2);
hold on;
a4=plot(X3(1,:), X3(2,:), '-','LineWidth', 2);
hold on;
a5=plot(Xa1(1,:), Xa1(2,:), '--','LineWidth', 2);
hold on;
a6=plot(Xa2(1,:), Xa2(2,:), '--','LineWidth', 2);
hold on;
a7=plot(Xa3(1,:), Xa3(2,:), '--','LineWidth', 2);
hold on;
% circle
thetaa=0:0.1:2*pi;
Circle1=X2(1,6)+0.3*cos(thetaa);
Circle2=X2(2,6)+0.3*sin(thetaa);
plot(Circle1,Circle2,'k--','linewidth',1);
Circle1=X3(1,6)+0.3*cos(thetaa);
Circle2=X3(2,6)+0.3*sin(thetaa);
plot(Circle1,Circle2,'k--','linewidth',1);
% marker
plot(X1(1,1), X1(2,1),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X2(1,1), X2(2,1),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X3(1,1), X3(2,1),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);

plot(X1(1,6), X1(2,6),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X2(1,6), X2(2,6),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X3(1,6), X3(2,6),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);

plot(X1(1,80), X1(2,80),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X2(1,80), X2(2,80),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X3(1,80), X3(2,80),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X1(1,151), X1(2,151),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X2(1,151), X2(2,151),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(X3(1,151), X3(2,151),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);

% marker
plot(Xa1(1,1), Xa1(2,1),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(Xa2(1,1), Xa2(2,1),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(Xa3(1,1), Xa3(2,1),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);

plot(Xa1(1,80), Xa1(2,80),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(Xa2(1,80), Xa2(2,80),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(Xa3(1,80), Xa3(2,80),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(Xa1(1,151), Xa1(2,151),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(Xa2(1,151), Xa2(2,151),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
plot(Xa3(1,151), Xa3(2,151),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
xlabel('x [m]','FontSize',12);
ylabel('y [m]','FontSize',12);
plot([X1(1,1), X2(1,1)],[X1(2,1), X2(2,1)],'r-.');
plot([X2(1,1), X3(1,1)],[X2(2,1), X3(2,1)],'r-.');
plot([X3(1,1), X1(1,1)],[X3(2,1), X1(2,1)],'r-.');

plot([X1(1,80), X2(1,80)],[X1(2,80), X2(2,80)],'r-.');
plot([X2(1,80), X3(1,80)],[X2(2,80), X3(2,80)],'r-.');
plot([X3(1,80), X1(1,80)],[X3(2,80), X1(2,80)],'r-.');

plot([X1(1,151), X2(1,151)],[X1(2,151), X2(2,151)],'r-.');
plot([X2(1,151), X3(1,151)],[X2(2,151), X3(2,151)],'r-.');
plot([X3(1,151), X1(1,151)],[X3(2,151), X1(2,151)],'r-.');

plot([Xa1(1,1), Xa2(1,1)],[Xa1(2,1), Xa2(2,1)],'b--');
plot([Xa2(1,1), Xa3(1,1)],[Xa2(2,1), Xa3(2,1)],'b--');
plot([Xa3(1,1), Xa1(1,1)],[Xa3(2,1), Xa1(2,1)],'b--');

plot([Xa1(1,80), Xa2(1,80)],[Xa1(2,80), Xa2(2,80)],'b--');
plot([Xa2(1,80), Xa3(1,80)],[Xa2(2,80), Xa3(2,80)],'b--');
plot([Xa3(1,80), Xa1(1,80)],[Xa3(2,80), Xa1(2,80)],'b--');

plot([Xa1(1,151), Xa2(1,151)],[Xa1(2,151), Xa2(2,151)],'b--');
plot([Xa2(1,151), Xa3(1,151)],[Xa2(2,151), Xa3(2,151)],'b--');
plot([Xa3(1,151), Xa1(1,151)],[Xa3(2,151), Xa1(2,151)],'b--');
grid on;
hold off;
legend([a1,a2,a3,a4,a5,a6,a7],{'reference','AUV1', 'AUV2', 'AUV3','AUVa1', 'AUVa2', 'AUVa3'});




figure(2)
subplot(3,1,1)
plot(TIme(1:Tstep),Ua1(1,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U1(1,1:Tstep),'LineWidth', 2);
hold on;
plot(TIme(1:Tstep),Ua2(1,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U2(1,1:Tstep),'LineWidth', 2);
hold on;
plot(TIme(1:Tstep),Ua3(1,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U3(1,1:Tstep),'LineWidth', 2);
grid on;
ylabel('u1 [N]','FontSize',12);

subplot(3,1,2)
plot(TIme(1:Tstep),Ua1(2,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U1(2,1:Tstep),'LineWidth', 2);
hold on;
plot(TIme(1:Tstep),Ua2(2,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U2(2,1:Tstep),'LineWidth', 2);
hold on;
plot(TIme(1:Tstep),Ua3(2,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U3(2,1:Tstep),'LineWidth', 2);
grid on;
ylabel('u2 [N]','FontSize',12);

subplot(3,1,3)
plot(TIme(1:Tstep),Ua1(3,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U1(3,1:Tstep),'LineWidth', 2);
hold on;
plot(TIme(1:Tstep),Ua2(3,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U2(3,1:Tstep),'LineWidth', 2);
hold on;
plot(TIme(1:Tstep),Ua3(3,1:Tstep), '--','LineWidth', 2);
hold on;
plot(TIme(1:Tstep),U3(3,1:Tstep),'LineWidth', 2);
hold on;
grid on;
legend({'AUVa1','AUV1','AUVa2','AUV2','AUVa3','AUV3'});
xlabel('Time [s]','FontSize',12);
ylabel('u3 [Nm]','FontSize',12);

rmpath('./utilis/','./ch6header/')