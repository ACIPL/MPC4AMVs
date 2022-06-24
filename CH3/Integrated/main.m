%% ------------------------------------------------------------------------
% License information: You are free to use or extend the code for research
% or education purposes provided that you provide clear a reference
% including:
% [1] Shen, C., Shi, Y. and Buckham, B., 2016. Integrated path planning and
%     tracking control of an AUV: A unified receding horizon optimization
%     approach. IEEE/ASME Transactions on Mechatronics.
% [2] Shi, Y., Shen, C., Wei, H., Zhang, K., 2022. Advanced model
%     predictive control for autonomous marine vehicels, Springer.
%
% Code description: Simulation of the Integrated path planning and MPC
% tracking control.
% Step 1) a feasible path is generated;
% Step 2) MPC generates the control input for the AUV to track the path.
% -------------------------------------------------------------------------

%% Offline setup 
addpath('./ch3header/','./utilis/')
clear;
clc;
load('SensorData.mat')

% AUV system parameters 
m=116; Iz=13.1;
X_udot=-167.6; Y_vdot=-477.2; N_rdot=-15.9;
Xu=26.9; Yv=35.8; Nr=3.5;
Du=241.3; Dv=503.8; Dr=76.9;

Mx=m-X_udot; My=m-Y_vdot; Mpsi=Iz-N_rdot;
coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr];

% Initialize MPC parameters
N = 8;
Tstep = 400;
Tstep2 = Tstep-N;

X0 = [0;1;0;0;0;0];
nx = length(X0);
nu = 3;

Xall = zeros(nx,Tstep+1);
Uall = zeros(nu,Tstep);
eta_Ref_all = zeros(3,Tstep);
Xall2 = zeros(nx,Tstep2+1);
Uall2 = zeros(nu,Tstep2);
Xall(:,1) = X0;
Xall2(:,1) = X0;
U0 = Uall2(:,1:N);

% the control input bound
Fu_max = 2000; Fv_max = 2000; Fr_max = 1000;
Umax0 = [Fu_max;Fv_max;Fr_max];
Umin0 = [-Fu_max;-Fv_max;-Fr_max];
for i = 1:1:N
    Umax(1+nu*(i-1):nu*i,1) = Umax0;
    Umin(1+nu*(i-1):nu*i,1) = Umin0;
end

% the weighting matrix of MPC
Q = 0.4 * eye(6); R = 0.01 * eye(3); Qf = 0.5 * eye(6);
weights = {Q, R, Qf};
K1 = 0.5 * eye(3); K2 = -1 * eye(3); K = [K1 K2];
disturbance=[0;0;0];
ndof = 3;
T = 0.1;

% initialize the controlled AUV
auv = AUV(coef, ndof, X0, U0);

% initialize the internal prediction model for MPC 
internal_auv = AUV(coef, ndof, X0, U0);
mpc1 = MPC_controller(N,internal_auv,weights,Umax0,Umin0);

% planning parameter setup
eps = 1e-1;
stroke = 20; % 20 * 10 cm = 2 m
y_scale = 0.1;

nsp = 40;
Nbrk = 3;
ord = 4;
mlti = 4;
SpAry = cell(nsp,1);
everdone = 0;
bbrks = everdone+0:1:everdone+Nbrk;
X = (everdone+0:1:everdone+Nbrk)';
Y = zeros(Nbrk+1,1);

for i = 1:1:Nbrk+1
    Y(i) = data1_together(1,100*(everdone+i-1)+1);
end
knots = augknt(bbrks,mlti);
nknt = length(knots);
nco = nknt-ord;

x0 = zeros(nco,1);
x0 = x0+1e-5;
plannar = Path_plannar(N,auv,stroke,mlti);

x = plannar.calc_path1(x0,bbrks,X,Y,T);
x = x';
sp = spmak(knots,x);

y0 = fnval(sp,everdone+1);
y0d = fnval(fnder(sp,1),everdone+1);
y02d = fnval(fnder(sp,2),everdone+1);
SpAry(1,1) = {sp};

iter_round = 1;
everdone = 1;

M = 10;
T_plt = T;
Brk_length = 1;
N_every_spline = Brk_length/T_plt;
N_total = nsp*N_every_spline;
% rescale the data into metric (in meter) | [xr;S(xr);S'(xr);S"(xr);S"'(xr)];
Spline_map2 = zeros(N_total,5);
[m,n] = size(Spline_map2);

% The upper and lower bound of operation space
X = (0:1:nsp)';
Y = X;
for i = 1:1:nsp+1
    Y(i) = data1_together(1,100*(i-1)+1);
end
x1 = X(1):T_plt:X(nsp+1);
x1 = x1';
nx1 = length(x1);
y1 = zeros(nx1,1);
y2 = y1;
for i = 1:1:nx1
    y1(i) = y_scale*plygfun1(x1(i),X,Y);
    y2(i) = y_scale*plygfun2(x1(i),X,Y,stroke);
end

u0 = zeros(nu*N,1);
% Xplus = X0;
P = zeros(nx+nu,N);
t = 0;
Xplus = X0;
%% ========================================================================
% ========================= Main function =================================
% =========================================================================
for i = 1:1:Tstep2
    % Planning
    if mod(i+N-1,M) == 0
        everdone = iter_round*1;
        if everdone == 30
            ord = 5;
            mlti = 5;
        end
        bbrks = everdone+0:1:everdone+Nbrk;
        X = (everdone+0:1:everdone+Nbrk)';
        
        for l = 1:1:Nbrk+1
            Y(l) = data1_together(1,100*(everdone+l-1)+1);
        end
        
        knots = augknt(bbrks,mlti);
        nknt = length(knots);
        nco = nknt-ord;
        
        x = plannar.calc_path2(x0,bbrks,X,Y,y0,y0d,y02d,eps,T);
        x = x';
        sp = spmak(knots,x);
        
        y0 = fnval(sp,everdone+1);
        y0d = fnval(fnder(sp,1),everdone+1);
        y02d = fnval(fnder(sp,2),everdone+1);
        SpAry(everdone+1,1) = {sp};
        iter_round = iter_round+1;
    end
    
    % Generate the path reference
    [eta_Ref_all(:,i),P] = plannar.reference(t,Spline_map2,SpAry,i,M,T_plt,y_scale,T);
    ord = 4;
    mlti = 4;
    
    % Calculate the control input for the AUV
    tic
    u = mpc1.calc_control(u0,X0,P,T);
    toc
    u_actual = u(1:nu,1);
    Uall2(:,i) = u_actual;
    
    % Update the system state
    auv.advance( Uall2(:,i), disturbance, T);
    Xplus = auv.X;
    Xall2(:,i+1) = Xplus;
    X0 = Xplus;
    
    % Construction of a feasible control u0
    Xplus2 = Xplus;
    
    for k = 1:1:N  % partition of u
        for j = 1:1:nu
            U(j,k) = u((k-1)*nu+j,1);
        end
    end
    
    for j=1:1:N-1 % state prediction
        Xplus2=internal_auv.dynamics_discrete( Xplus2,U(:,j),T );
    end
    
    XN = Xplus2;
    XrN = P(1:6,N);
    XeN = internal_auv.ErrorState( XN, XrN );
    u0=[u(nu+1:nu*N,1);K*XeN];
    
    t=t+T;
end

%% ========================================================================
% ============================ Plot =======================================
% =========================================================================

% plot the trajectory
figure(1)
plot(x1(:,1),y1(:,1),'g','LineWidth',2);
hold on;
plot(x1(:,1),y2(:,1),'g','LineWidth',2);
grid on;
plot(eta_Ref_all(1,1:Tstep2),eta_Ref_all(2,1:Tstep2),'k','LineWidth',2);
hold on;
plot(Xall2(1,1:Tstep2),Xall2(2,1:Tstep2),'r','LineWidth',2);
hold on;

% plot the control inputs
figure(2)
subplot(3,1,1)
plot(Uall2(1,:),'r');
hold on;
grid on;
subplot(3,1,2)
plot(Uall2(2,:),'r');
hold on;
grid on;
subplot(3,1,3)
plot(Uall2(3,:),'r');
hold on;
grid on;

rmpath('./ch3header/','./utilis/')