%% ------------------------------------------------------------------------
% License information: You are free to use or extend the code for research
% or education purposes provided that you provide clear a reference
% including:
% [1] Wei, Henglai, et al. "Robust distributed model predictive platooning 
%     control for heterogeneous autonomous surface vehicles." Control 
%     Engineering Practice 107 (2021): 104655.
% [2] Shi, Y., Shen, C., Wei, H., Zhang, K., 2022. Advanced model
%     predictive control for autonomous marine vehicels, Springer.
%
% Code description: 
% Simulation of DMPC for ASVs
%% ------------------------------------------------------------------------

addpath('./utilis/','./ch7header/')
clc;
clear all;

%======================================================================%
% System coefficients
%======================================================================%

m0 = 15; X_udot0 = -2; Xu0 = 0.72; Du0 = 2.05; Mx0 = m0-X_udot0; 
m1 = 18; X_udot1 = -2; Xu1 = 0.8; Du1 = 2.0; Mx1 = m1-X_udot1;
m2 = 12; X_udot2 = -3; Xu2 = 1.3; Du2 = 1.6; Mx2 = m2-X_udot2;
m3 = 14; X_udot3 = -1.6; Xu3 = 1.7; Du3 = 2.3; Mx3 = m3-X_udot3;
m4 = 13.7; X_udot4 = -2.6; Xu4 = 1.35; Du4 = 1.8; Mx4 = m4-X_udot4;
m5 = 13; X_udot5 = -1.3; Xu5 = 1.0; Du5 = 1.3; Mx5 = m5-X_udot5;
m6 = 14; X_udot6 = -3.6; Xu6 = 1.2; Du6 = 1.3; Mx6 = m6-X_udot6;
disturbance = 0;

%======================================================================%
% create 3 ASVs 
%======================================================================%
coef0 = [m0; X_udot0; Xu0; Du0];
coef1 = [m1; X_udot1; Xu1; Du1];
coef2 = [m2; X_udot2; Xu2; Du2];
coef3 = [m3; X_udot3; Xu3; Du3];
coef4 = [m4; X_udot4; Xu4; Du4];
coef5 = [m5; X_udot5; Xu5; Du5];
coef6 = [m6; X_udot6; Xu6; Du6];
ndof = 1;
X0_0 = [0; 1];
X1_0 = [-2.5; 1];
X2_0 = [-5; 1];
X3_0 = [-7.5; 1];
X4_0 = [-10; 1];
X5_0 = [-12.5; 1];
X6_0 = [-15; 1];
U0 = 0; 
% Denote the actual system
asv0 = ASV(coef0, ndof, X0_0, U0);
asv1 = ASV(coef1, ndof, X1_0, U0);
asv2 = ASV(coef2, ndof, X2_0, U0);
asv3 = ASV(coef3, ndof, X3_0, U0);
asv4 = ASV(coef4, ndof, X4_0, U0);
asv5 = ASV(coef5, ndof, X5_0, U0);
asv6 = ASV(coef6, ndof, X6_0, U0);
% Denote the prediction model
internal_asv0 = ASV(coef0, ndof, X0_0, U0);
internal_asv1 = ASV(coef1, ndof, X1_0, U0);
internal_asv2 = ASV(coef2, ndof, X2_0, U0);
internal_asv3 = ASV(coef3, ndof, X3_0, U0);
internal_asv4 = ASV(coef4, ndof, X4_0, U0);
internal_asv5 = ASV(coef5, ndof, X5_0, U0);
internal_asv6 = ASV(coef6, ndof, X6_0, U0);

%======================================================================%
% create three sine2D reference trajectories
%======================================================================%
% the initial velicity, acceleration and time instant
v0 = 1; a0 = 0; t0 = 0; 
% ref_trajectory = line1D(v0, a0, t0); 
% ref = [ref_trajectory.xR; ref_trajectory.xRd; ref_trajectory.xRdd];
d0r = [0; 0]; d1r = [-2.5; 0]; d2r = [-5; 0];
d3r = [-7.5; 0]; d4r = [-10; 0]; d5r = [-12.5; 0]; d6r = [-15; 0];

%==========================================================================
% Other design parameters
%==========================================================================
% sampling time
dt = 0.1;
% simulation step
Tstep = 60;
nx = length(X1_0);
nu = length(U0);

%======================================================================%
% create Robust DMPC platooning controllers
%======================================================================%
% Some intial setup for RDMPC
N = 6;  % prediction horizon N
M = nx + nu;
t = 0;
[u0_0, u1_0, u2_0, u3_0, u4_0, u5_0, u6_0] = deal(zeros(nu*N,1));
% Weighting Matrices
Qi = diag([1e3 1e2]);
Ri = 1e-4 ;
Qfi = diag([1e4 1e3]);
Fi = diag([1e2 1e2]);
Hi = diag([1e2 1e1]);
weights = {Qi,Ri,Qfi,Fi,Hi};

U0_max = 20.67;
U1_max = 22.16;
U2_max = 21.32;
U3_max = 22.93;
U4_max = 26.24;
U5_max = 26.02;
U6_max = 21.38;

rdmpc0 = RDMPC_controller(N,internal_asv0,weights,U0_max,-U0_max);
rdmpc1 = RDMPC_controller(N,internal_asv1,weights,U1_max,-U1_max);
rdmpc2 = RDMPC_controller(N,internal_asv2,weights,U2_max,-U2_max);
rdmpc3 = RDMPC_controller(N,internal_asv3,weights,U3_max,-U3_max);
rdmpc4 = RDMPC_controller(N,internal_asv4,weights,U4_max,-U4_max);
rdmpc5 = RDMPC_controller(N,internal_asv5,weights,U5_max,-U5_max);
rdmpc6 = RDMPC_controller(N,internal_asv6,weights,U6_max,-U6_max);
%======================================================================%
% DLMPC simulation
%======================================================================%
% To store the simulation data
[X0, X1, X2, X3, X4, X5, X6, Ref] = deal(zeros(nx,Tstep+1));
[U0, U1, U2, U3, U4, U5, U6] = deal(zeros(nu,Tstep));

[Xreal0, Xreal1, Xreal2, Xreal3, Xreal4, Xreal5, Xreal6] = deal(zeros(nx,Tstep+1));
Ref0 = [0;1];
X0(:,1) = X0_0; X1(:,1) = X1_0; X2(:,1) = X2_0; X3(:,1) = X3_0;
X4(:,1) = X4_0; X5(:,1) = X5_0; X6(:,1) = X6_0; Ref(:,1) = Ref0;
Xreal0(:,1) = X0_0; Xreal1(:,1) = X1_0; Xreal2(:,1) = X2_0; 
Xreal3(:,1) = X3_0; Xreal4(:,1) = X4_0; Xreal5(:,1) = X5_0; 
Xreal6(:,1) = X6_0; 
% initial broadcast information
[X0_b, X1_b, X2_b, X3_b, X4_b, X5_b, X6_b] = deal(zeros(nx,N));
% to store the time sequence
TIme = zeros(Tstep+1,1); 
tic
for i=1:1:Tstep
    tic
    t1 = t;
    TIme(i) = t1;
    %======================================================================
    % ASV 0
    %======================================================================
    P = RefGen(nx, N, t1, dt, Ref0);
    Ref0 = P(:,1);
    Ref(:,i) = P(:,1);
    % Reference for ASV 0
    P0 = P + d0r;
    X0_neighbors = {X0_b};
    [u0, X0_b] = rdmpc0.calc_control(P0, X0_neighbors, X0_0, u0_0, dt, d0r);
    u0_actual = u0(1);
    U0(:,i) = u0_actual;
    % generate the norminal state
    X0(:,i+1) = X0_b(:,1);
    X0_0 = X0_b(:,1);
    % general the actual system state
    disturbance0 = 5*sin(0.7*pi*t);
    Xreal0(:,i+1) = asv0.dynamics_discrete(X0_0,u0_actual+disturbance0,dt);
    % Construction of a feasible control u0_0
    u0_0 = [u0(nu+1:nu*N,1);u0(nu*(N-1)+1:nu*N,1)];  
    
    %======================================================================
    % ASV 1
    %======================================================================
    % Reference for ASV 1
    P1 = P + d1r;
    X1_neighbors = {X0_b, X1_b};
    [u1, X1_b] = rdmpc1.calc_control(P1, X1_neighbors, X1_0, u1_0, dt, d1r);
    
    u1_actual = u1(1:nu, 1); 
    U1(:,i) = u1_actual;
    % generate the norminal state
    X1(:,i+1) = X1_b(:,1)+d1r;
    X1_0 = X1_b(:,1)+d1r;
    % general the actual system state
    disturbance1 =  4.5*sin(0.5*pi*t+pi/2);
    Xreal1(:,i+1) = asv1.dynamics_discrete(X1_0,u1_actual+disturbance1,dt);
    % Construction of a feasible control u1_0
    u1_0 = [u1(nu+1:nu*N,1);u1(nu*(N-1)+1:nu*N,1)];  
    
    %======================================================================
    % ASV 2 
    %======================================================================
    P2 = P + d2r;
    X2_neighbors = {X0_b, X1_b};
    [u2, X2_b] = rdmpc2.calc_control(P2, X2_neighbors, X2_0, u2_0, dt, d2r);
    % The control input apply to the ASV and store it 
    u2_actual = u2(1:nu, 1); 
    U2(:,i) = u2_actual;
    % The norminal state
    X2(:,i+1) = X2_b(:,1)+d2r;
    X2_0 = X2_b(:,1)+d2r;
    % general the actual system state
    disturbance2 = 5.5*sin(0.8*pi*t+pi/3);
    Xreal2(:,i+1) = asv2.dynamics_discrete(X2_0,u2_actual+disturbance2,dt);
    %======================================================================%
    % Construction of a feasible control u2_0 
    u2_0 = [u2(nu+1:nu*N,1);u2(nu*(N-1)+1:nu*N,1)]; 
    
    %======================================================================
    % AUV 3
    %======================================================================
    P3 = P + d3r;
    X3_neighbors = {X0_b, X2_b};
    [u3, X3_b] = rdmpc3.calc_control(P3, X3_neighbors, X3_0, u3_0, dt, d3r);
    % The control input apply to the ASV 3 and store it 
    u3_actual = u3(1:nu, 1); 
    U3(:,i) = u3_actual;
    % The norminal state
    X3(:,i+1) = X3_b(:,1)+d3r;
    X3_0 = X3_b(:,1)+d3r;
    % general the actual system state
    disturbance3 =  5.3*sin(0.75*pi*t-pi/4);
    Xreal3(:,i+1) = asv3.dynamics_discrete(X3_0,u3_actual+disturbance3,dt);
    % Construction of a feasible control u1_0 
    u3_0 = [u3(nu+1:nu*N,1);u3(nu*(N-1)+1:nu*N,1)]; 
    
    %======================================================================
    % ASV 4
    %======================================================================
    P4 = P + d4r;
    X4_neighbors = {X0_b, X3_b};
    [u4, X4_b] = rdmpc4.calc_control(P4, X4_neighbors, X4_0, u4_0, dt, d4r);
    % The control input apply to the ASV and store it 
    u4_actual = u4(1:nu,1); 
    U4(:,i) = u4_actual;
    % The norminal state
    X4(:,i+1) = X4_b(:,1)+d4r;
    X4_0 = X4_b(:,1)+d4r;
    % general the actual system state
    disturbance4 =  6.1*sin(0.9*pi*t+pi/5);
    Xreal4(:,i+1) = asv4.dynamics_discrete(X4_0,u4_actual+disturbance4,dt);
    % Construction of a feasible control u1_0 
    u4_0 = [u4(nu+1:nu*N,1);u4(nu*(N-1)+1:nu*N,1)]; 
    
    %======================================================================
    % ASV 5
    %======================================================================
    P5 = P + d5r;
    X5_neighbors = {X0_b, X4_b};
    [u5, X5_b] = rdmpc5.calc_control(P5, X5_neighbors, X5_0, u5_0, dt, d5r);
    % Store the calculated control input
    u5_actual = u5(1:nu, 1); 
    U5(:,i) = u5_actual;
    % Store the norminal state
    X5(:,i+1) = X5_b(:,1)+d5r;
    X5_0 = X5_b(:,1)+d5r;
    % general the actual system state
    disturbance5 =  5.89*sin(0.8*pi*t+2*pi/3);
    Xreal5(:,i+1) = asv5.dynamics_discrete(X5_0,u5_actual+disturbance5,dt);
    %======================================================================%
    % Construction of a feasible control u1_0 
    u5_0 = [u5(nu+1:nu*N,1);u5(nu*(N-1)+1:nu*N,1)]; 
    
    %======================================================================
    % ASV 6
    %======================================================================
    P6 = P + d6r;
    X6_neighbors = {X0_b, X5_b};
    [u6, X6_b] = rdmpc6.calc_control(P6, X6_neighbors, X6_0, u6_0, dt, d6r);
    % Store the calculated control input
    u6_actual = u6(1:nu, 1); 
    U6(:,i) = u6_actual;
    % Store the norminal state
    X6(:,i+1) = X6_b(:,1)+d6r;
    X6_0 = X6_b(:,1)+d6r;
    % general the actual system state
    disturbance6 = 6*sin(0.65*pi*t-pi/6);
    Xreal6(:,i+1) = asv1.dynamics_discrete(X6_0,u6_actual+disturbance6,dt);
    %======================================================================%
    % Construction of a feasible control u1_0 
    u6_0 = [u6(nu+1:nu*N,1);u6(nu*(N-1)+1:nu*N,1)]; 
    i
    t = t + dt;
    toc
end

%==========================================================================
plot

rmpath('./utilis/','./ch7header/')