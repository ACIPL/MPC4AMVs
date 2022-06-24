%% ------------------------------------------------------------------------
% License information: You are free to use or extend the code for research
% or education purposes provided that you provide clear a reference
% including:
% [1] Wei, Henglai, Chao Shen, and Yang Shi. "Distributed Lyapunov-based
%     model predictive formation tracking control for autonomous underwater
%     vehicles subject to disturbances." IEEE Transactions on Systems, Man,
%     and Cybernetics: Systems 51.8 (2021): 5198-5208.
% [2] Shi, Yang, Shen, Chao, Wei, Henglai, and Zhang, Kunwu, 2022. Advanced model
%     predictive control for autonomous marine vehicels, Springer.
%
% Code description: Test 3, DLMPC with ESO: tracking sin shape

clear;
close all;
clc;
addpath('./utilis/')
global rho R  Qi Ri Pi beta
Qi = diag([10^4 10^4 10^3 10^2 100 100]);
Ri = diag([10^(-3) 10^-3 10^-3]);
Pi = diag([10^4 10^4 10^3 10^2 10^2 10^2]);
rho = 1;
R = diag([rho rho]);
beta = 10^4;

% parameter for formation and tracking
global d1r d2r d3r d12 d21 d13 d31 d23 d32
d1r = [0; 1; 0]; d2r = [0; -1; 0]; d3r = [-1; 0; 0];
d12 = [0, 2]; d21 = [0, -2]; d13 = [1, 1]; d31 = [-1, -1];
d23 = [1, -1]; d32 = [-1, 1];

mpciterations = 150;
N             = 5;
T             = 0.1;
tol_opt       = 1e-8;
opt_option    = 0;
iprint        = 0;
type          = 'differential equation';
atol_ode_real = 1e-12;
rtol_ode_real = 1e-12;
atol_ode_sim  = 1e-4;
rtol_ode_sim  = 1e-4;
iprint = 0;

t = 0; mpciter = 0;
tmeasure1 = 0.0; tmeasure2 = 0.0; tmeasure3 = 0.0;
xmeasure1 = [0; 2.5; pi/2; 0.3; 0; 0];
xmeasure2 = [-0.5; -0.5; 3*pi/4; 0.4; 0; 0];
xmeasure3 = [-1.5; 1.1; -pi/2; 0.3; 0; 0];
u10 = 0.01*ones(3,N); u20= 0.01*ones(3,N); u30 = 0.01*ones(3,N);

t1 = tmeasure1; t2 = tmeasure2; t3 = tmeasure3;
x1 = xmeasure1; x2 = xmeasure2; x3 = xmeasure3;
u1 = []; u2 = []; u3 = [];
X1 = []; X2= []; X3 = [];
X1 = x1'; X2 = x2'; X3 = x3';
xhat1 = zeros(N+1,6); xhat2 = zeros(N+1,6); xhat3 = zeros(N+1,6);
hat_dxa1 = zeros(9,1);hat_dxa2 = zeros(9,1);hat_dxa3 = zeros(9,1);
Ua1 = [];Ua2 = [];Ua3 = [];
xa1 = xmeasure1; xa2 = xmeasure2; xa3 = xmeasure3;
Xa1 = [];Xa2 = [];Xa3 = [];
Xa1 = xa1'; Xa2 = xa2'; Xa3 = xa3'; % states under auxiliary law

hat_xa1(1:6) = xa1; hat_xa1(7:9) = rand(3,1); hat_xa1 = hat_xa1';
hat_xa2(1:6) = xa2; hat_xa2(7:9) = rand(3,1); hat_xa2 = hat_xa2';
hat_xa3(1:6) = xa3; hat_xa3(7:9) = rand(3,1); hat_xa3 = hat_xa3';
hat_Xa1 = []; hat_Xa2 = []; hat_Xa3 = []; % estimated states under auxiliary law
hat_Xa1 = hat_xa1'; hat_Xa2 = hat_xa2'; hat_Xa3 = hat_xa3';
lv1 = []; lv2=[]; lv3=[]; % Lyapunov function value under auxiliary law
mv1= []; mv2=[]; mv3=[]; % Lyapunov function value under MPC
hat_W1 = []; hat_W2 = []; hat_W3 = []; % estimated disturbances
D1=[];D2=[];D3=[]; % real disturbances
E1 = []; E2 = []; E3 = []; % disturbance error
m=116;
Iz=13.1;
X_udot=-167.6;
Y_vdot=-477.2;
N_rdot=-15.9;
Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;
M=diag([Mx;My;Mpsi]);
while (mpciter < mpciterations)
    % ESO-based auxiliary control law for AUV1
    [ua1,V1,LV1] = bs2(t, hat_xa1,d1r);
    Ua1 = [Ua1; ua1'];
    % to estimate the states
    [hat_dxa1,dxa1,hat_w1,d1,ee1] = ESO(t,hat_xa1,xa1,ua1);
    xa11 = hat_xa1;
    xa1 = xa1 + dxa1*T;
    Xa1 = [Xa1; xa1'];
    hat_xa1 = hat_xa1 + hat_dxa1*T;
    hat_Xa1 = [hat_Xa1; hat_xa1'];
    hat_xa1(7:9)=hat_w1;
    hat_W1 = [hat_W1;hat_w1'];
    lv1 = [lv1;LV1];
    D1 = [D1;d1'];
    E1 = [E1;ee1];
    
    
    % DLMPC for AUV1
    [tmeasure1, xmeasure1, u, xhat10]=ddmpc(@runningcosts1, @terminalcosts1,...
        @constraints1,@terminalconstraints1, @linearconstraints1, @system_ct, ...
        N, T, t, x1', u10, xhat2, xhat3, V1, tol_opt, opt_option, ...
        type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim,iprint);
    [~,dx1,~,~,~] = ESO(t,xa11,x1,u(:,1));
    [V11,MV1] = bs22(t,x1',u(:,1),d1r);
    x1 = x1 + dx1*T;
    t1 = [t1;tmeasure1];
    X1 = [X1;x1'];
    u1 = [u1;u'];
    mv1 = [mv1;MV1];
    
    
    %ESO-based auxiliary control law for AUV2
    [ua2,V2,LV2] = bs2(t, hat_xa2,d2r);
    Ua2 = [Ua2; ua2']; % control inputs
    [hat_dxa2,dxa2,hat_w2,d2,ee2] = ESO(t,hat_xa2,xa2,ua2);
    xa22 = hat_xa2;
    xa2 = xa2 + dxa2*T;
    Xa2 = [Xa2; xa2']; % real states
    hat_xa2 = hat_xa2 + hat_dxa2*T;
    hat_Xa2 = [hat_Xa2; hat_xa2']; %estimated states
    hat_xa2(7:9)=hat_w2;
    hat_W2 = [hat_W2;hat_w2'];
    lv2 = [lv2;LV2];
    D2 = [D2;d2'];
    E2 = [E2;ee2];
    
    
    %DLMPC for AUV2
    [tmeasure2, xmeasure2, u, xhat20]=ddmpc(@runningcosts2, @terminalcosts2,...
        @constraints2, @terminalconstraints2, @linearconstraints2, @system_ct, ...
        N, T, t, x2', u20, xhat1, xhat3, V2, tol_opt, opt_option, ...
        type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim,iprint);
    [~,dx2,~,~,~] = ESO(t,xa22,x2,u(:,1));
    [V22,MV2] = bs22(t,x2',u(:,1),d2r);
    x2 = x2 + dx2*T;
    X2 = [X2;x2'];
    u2 = [u2;u'];
    mv2=[mv2;MV2];
    
    
    %ESO-based auxiliary control law for AUV3
    [ua3,V3,LV3] = bs2(t,hat_xa3,d3r);
    Ua3 = [Ua3; ua3'];
    [hat_dxa3,dxa3,hat_w3,d3,ee3] = ESO(t,hat_xa3,xa3,ua3);
    xa33 = hat_xa3;
    hat_xa3 = hat_xa3 + hat_dxa3*T;
    hat_Xa3 = [hat_Xa3; hat_xa3'];
    xa3 = xa3 + dxa3*T;
    hat_xa3(7:9)=hat_w3;
    hat_W3 = [hat_W3;hat_w3'];
    Xa3 = [Xa3; xa3'];
    lv3 = [lv3;LV3];
    D3 = [D3;d3'];
    E3 = [E3;ee3];
    
    
    %DLMPC for AUV3
    [tmeasure3, xmeasure3, u, xhat30]=ddmpc(@runningcosts3, @terminalcosts3,...
        @constraints3, @terminalconstraints3, @linearconstraints3, @system_ct, ...
        N, T, t, x3', u30, xhat1, xhat2, V3, tol_opt, opt_option, ...
        type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim,iprint);
    [~,dx3,~,~,~] = ESO(t,xa33,x3,u(:,1));
    [V33,MV3] = bs22(t,x3',u(:,1),d3r);
    x3 = x3 + dx3*T;
    X3 = [X3;x3'];
    u3 = [u3;u'];
    mv3=[mv3;MV3];
    
    xhat1=xhat10;
    xhat2=xhat20;
    xhat3=xhat30;
    
    mpciter = mpciter+1
    t = t+T;
end

rmpath('./utilis/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUV1 %yd+1
%%%%%%%%%%%%%%%%%%%% Definition of the cost functions %%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts1(t, x, u, xhat11, xhat22)
global Qi Ri d12 d13 beta
xd = 0.5*t;
yd = sin(0.5*t+pi/2);
psid = atan2(cos(0.5*t+pi/2),1);
ud = sqrt((0.5^2)+(0.5*cos(0.5*t+pi/2))^2);
vd = 0;
rd = -0.125*sin(0.5*t+pi/2)/ud^2;
xe = x-[xd yd+1 psid ud vd rd];

cost = xe * Qi * xe' + u' * Ri * u + ...
    (x(1,1:2) - xhat11(1,1:2) - d12) * diag([beta,beta])*...
    (x(1,1:2) - xhat11(1,1:2) - d12)' + (x(1,1:2) - xhat22(1,1:2)...
    -d13) * diag([beta,beta]) * (x(1,1:2) - xhat22(1,1:2) - d13)';
end

function cost = terminalcosts1(t, x, xhat11, xhat22)
global Pi d12 d13 beta
xd = 0.5*t;
yd = sin(0.5*t+pi/2);
psid = atan2(cos(0.5*t+pi/2),1);
ud = sqrt((0.5^2)+(0.5*cos(0.5*t))^2);
vd = 0;
rd = -0.125*sin(0.5*t+pi/2)/ud^2;
xe = x-[xd yd+1 psid ud vd rd];
cost = xe * Pi * xe' +...
    (x(1,1:2) - xhat11(1,1:2) - d12) * diag([beta,beta])*...
    (x(1,1:2) - xhat11(1,1:2) - d12)' + (x(1,1:2) - xhat22(1,1:2)...
    -d13) * diag([beta,beta]) * (x(1,1:2) - xhat22(1,1:2) - d13)';
end
%%%%%%%%%%%%%%%%%%% Definition of the cost constraints %%%%%%%%%%%%%%%%%%%%
function [c,ceq] = constraints1(t, x, u, V1)
global d1r
c = [];
ceq  = [];
[V11,MV11] = bs22(t, x, u,d1r);
c = V11 - 1.1*V1;
end
function [A, b, Aeq, beq, lb, ub] = linearconstraints1(t, x, ~)
A   = [];
b   = [];
Aeq = [];
beq = [];
lb = [-1000 -1000 -1000];
ub = [1000 1000 1000];
end

function [c,ceq] = terminalconstraints1(t, x)
c = [];
ceq  = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AUV2
%%%%%%%%%%%%%%%%%%%% Definition of the cost functions %%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts2(t, x, u, xhat11, xhat22)
global Qi Ri d23 d21 beta
xd = 0.5*t;
yd = sin(0.5*t+pi/2);
psid = atan2(cos(0.5*t+pi/2),1);
ud = sqrt((0.5^2)+(0.5*cos(0.5*t+pi/2))^2);
vd = 0;
rd = -0.125*sin(0.5*t+pi/2)/ud^2;
xe = x-[xd yd-1 psid ud vd rd];
cost = xe * Qi * xe' + u' * Ri * u + ...
    (x(1,1:2) - xhat11(1,1:2) - d21) * diag([beta,beta])*...
    (x(1,1:2) - xhat11(1,1:2) - d21)';
end

function cost = terminalcosts2(t, x, xhat11, xhat22)
global Pi d23 d21 beta
xd = 0.5*t;
yd = sin(0.5*t+pi/2);
psid = atan2(cos(0.5*t+pi/2),1);
ud = sqrt((0.5^2)+(0.5*cos(0.5*t+pi/2))^2);
vd = 0;
rd = -0.125*sin(0.5*t+pi/2)/ud^2;
xe = x-[xd yd-1 psid ud vd rd];
cost = xe * Pi * xe'+ ...
    (x(1,1:2) - xhat11(1,1:2) - d21) * diag([beta,beta])*...
    (x(1,1:2) - xhat11(1,1:2) - d21)';
end

function [c,ceq] = constraints2(t, x, u,V2)
global d2r
c = [];
ceq  = [];
[V22,MV22] = bs22(t, x, u,d2r);
c = V22 - 1.1*V2;
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints2(t, x,~)
A   = [];
b   = [];
Aeq = [];
beq = [];
lb = [-1000 -1000 -1000];
ub = [1000 1000 1000];
end

function [c,ceq] = terminalconstraints2(t, x)
c = [];
ceq  = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUV3
%%%%%%%%%%%%%%%%%%%% Definition of the cost functions %%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts3(t, x, u, xhat11, xhat22)
global Qi Ri d32 d31 beta
xd = 0.5*t;
yd = sin(0.5*t+pi/2);
psid = atan2(cos(0.5*t+pi/2),1);
ud = sqrt((0.5^2)+(0.5*cos(0.5*t+pi/2))^2);
vd = 0;
rd = -0.125*sin(0.5*t+pi/2)/ud^2;
xe = x-[xd-1 yd psid ud vd rd];
cost = xe * Qi * xe' + u' * Ri * u + ...
    (x(1,1:2) - xhat11(1,1:2) - d31) * diag([beta,beta])*...
    (x(1,1:2) - xhat11(1,1:2) - d31)' ;
end

function cost = terminalcosts3(t, x, xhat11, xhat22)
global Pi d32 d31 beta
xd = 0.5*t;
yd = sin(0.5*t+pi/2);
psid = atan2(cos(0.5*t+pi/2),1);
ud = sqrt((0.5^2)+(0.5*cos(0.5*t+pi/2))^2);
vd = 0;
rd = -0.125*sin(0.5*t+pi/2)/ud^2;
xe = x-[xd-1 yd psid ud vd rd];
cost = xe * Pi * xe'+ ...
    (x(1,1:2) - xhat11(1,1:2) - d31) * diag([beta,beta])*...
    (x(1,1:2) - xhat11(1,1:2) - d31)';
end

function [c,ceq] = constraints3(t, x, u,V3)
global d3r
c =[];
ceq  = [];
[V33,MV33] = bs22(t, x, u,d3r);
c = V33 - 1.1*V3;

end

function [A, b, Aeq, beq, lb, ub] = linearconstraints3(t, x,~)
A   = [];
b   = [];
Aeq = [];
beq = [];
lb = [-1000 -1000 -1000];
ub = [1000 1000 1000];
end

function [c,ceq] = terminalconstraints3(t, x)
c = [];
ceq  = [];
end

% AUV model
function dx = system_ct(t,x,u,T)
% System coefficients
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
%
dx  = zeros(6,1);
dx(1) = x(4)*cos(x(3)) - x(5)*sin(x(3));
dx(2) = x(4)*sin(x(3)) + x(5)*cos(x(3));
dx(3) = x(6);
dx(4) = (My/Mx)*x(5)*x(6)-(Xu/Mx)*x(4)-(Du/Mx)*x(4)*abs(x(4))+u(1)/Mx;
dx(5) = -(Mx/My)*x(4)*x(6)-(Yv/My)*x(5)-(Dv/My)*x(5)*abs(x(5))+u(2)/My;
dx(6) = ((Mx-My)/Mpsi)*x(4)*x(5)-(Nr/Mpsi)*x(6)-(Dr/Mpsi)*x(6)*abs(x(6))+u(3)/Mpsi;
end

function [hat_dx,dx,hat_w,d,ee] = ESO(t,hat_x,x,u)
% Input: Time instant, estimated state, real state, control input;
% Output: real state (unknown), estimated state, disturbances;
% System coefficients for AUV
m=116;
Iz=13.1;
X_udot=-167.6;
Y_vdot=-477.2;
N_rdot=-15.9;
Xu=26.9;
Yv=35.8;
Nr=3.5;
Du=241.3;
Dv=503.8;
Dr=76.9;
Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;
M=diag([Mx;My;Mpsi]);
% updatea the control input
d = zeros(3,1);
d(1) = (0.2*x(5)*x(5)*x(5)+0.2*sin(0.7*t));
d(2) = (0.3*x(4)*x(6)+0.3*x(4)+0.1*sin(0.6*t));
d(3) = (-0.2*x(5)*x(5)-0.4*x(5)*x(4)-0.1*sin(0.9*t));
dx = zeros(6,1);
dx(1) = x(4)*cos(x(3)) - x(5)*sin(x(3));
dx(2) = x(4)*sin(x(3)) + x(5)*cos(x(3));
dx(3) = x(6);
dx(4) = (My/Mx)*x(5)*x(6)-(Xu/Mx)*x(4)-(Du/Mx)*x(4)*abs(x(4))+u(1)/Mx+d(1);
dx(5) = -(Mx/My)*x(4)*x(6)-(Yv/My)*x(5)-(Dv/My)*x(5)*abs(x(5))+(u(2))/My+d(2);
dx(6) = ((Mx-My)/Mpsi)*x(4)*x(5)-(Nr/Mpsi)*x(6)-(Dr/Mpsi)*x(6)*abs(x(6))+(u(3))/Mpsi+d(3);
% ESO
% system with disturbance
e1 = x(1:3) - hat_x(1:3);
e2 = x(4:6) - hat_x(4:6);
e3 = hat_x(7:9)-d;
ee = e1'*e1+e2'*e2+e3'*e3;
K1 = 10;
K2 = 10;
K3 = 20;

R2 = Rot(hat_x(3));
hat_dx  = zeros(9,1);
hat_dx(1) = hat_x(4)*cos(hat_x(3)) - hat_x(5)*sin(hat_x(3))+K1*e1(1);
hat_dx(2) = hat_x(4)*sin(hat_x(3)) + hat_x(5)*cos(hat_x(3))+K1*e1(2);
hat_dx(3) = hat_x(6)+K1*e1(3);

hat_dx(4) = (My/Mx)*x(5)*x(6)-(Xu/Mx)*x(4)-(Du/Mx)*x(4)*abs(x(4))+u(1)/Mx+hat_x(7);
hat_dx(5) = -(Mx/My)*x(4)*x(6)-(Yv/My)*x(5)-(Dv/My)*x(5)*abs(x(5))+u(2)/My+hat_x(8);
hat_dx(6) = ((Mx-My)/Mpsi)*x(4)*x(5)-(Nr/Mpsi)*x(6)-(Dr/Mpsi)*x(6)*abs(x(6))+u(3)/Mpsi+hat_x(9);
hat_dx(4:6) = hat_dx(4:6) + K2*e1;
hat_dx(7:9) = K3*e1;
hat_w = hat_dx(7:9);
end

function [ua2,V2,LV2] = bs2(t, x,dir)
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
Lambda=diag([1;1;1]);
PI=[0,-1,0;1,0,0;0,0,0];
Kp=diag([1;1;1]);
Kd=diag([1;1;1]);


xr = 0.5*t;
xrdot = 0.5;
xrddot = 0;
xrdddot = 0;

yr = sin(0.5*t+pi/2);
yrdot = 0.5*cos(0.5*t+pi/2);
yrddot = -0.5^2*sin(0.5*t+pi/2);
yrdddot = 0.5^3*cos(0.5*t+pi/2);

thetar = atan2(yrdot,xrdot);
thetardot = (xrdot*yrddot-yrdot*xrddot)/(xrdot^2+yrdot^2);
thetarddot=(1/(xrdot^2+yrdot^2)^2)*((xrddot*yrddot+xrdot*yrdddot-yrddot...
    *xrddot-yrdot*xrdddot)*(xrdot^2+yrdot^2)-(2*xrdot*xrddot+2*yrdot...
    *yrddot)*(xrdot*yrddot-yrdot*xrddot));

eta_r = [xr;yr;thetar]+dir;
eta_r_dot = [xrdot;yrdot;thetardot];
eta_r_ddot = [xrddot;yrddot;thetarddot];

% using backstep technique to contruct the control law
eta = x(1:3);
vel = x(4:6);
theta = x(3);
tilde_eta = eta - eta_r;
eta_d_dot = eta_r_dot - Lambda*tilde_eta;
R = Rot(theta);
eta_dot = R*vel;

u=vel(1);
v=vel(2);
r=vel(3);
C1=[0;0;My*v];
C2=[0;0;-Mx*u];
C3=[-My*v;Mx*u;0];
C=[C1,C2,C3];
D=diag([-Xu-Du*abs(u);-Yv-Dv*abs(v);-Nr-Dr*abs(r)]);
D=-D;

v_d = R'*eta_d_dot;
S = eta_dot - eta_d_dot;
Mstar= R*M*R';
Cstar=R*(C-M*R'*vel(3)*R*PI)*R';
Dstar=R*D*R';
thetadot = eta_dot(3);
v_d_dot=-thetadot*R'*PI*eta_d_dot+R'*(eta_r_ddot- Lambda*S);

tau=M*v_d_dot+C*v_d+D*v_d-R'*Kp*tilde_eta-R'*Kd*S-M*x(7:9);
V2=S'*R*(tau-M*v_d_dot-C*v_d+D*v_d+R'*Kp*tilde_eta-R'*Kd*S)-S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta;
LV2 = (5*tilde_eta'*Kp*tilde_eta+0.01*S'*Mstar*S);
ua2 = M*v_d_dot+C*v_d+D*v_d-R'*Kp*tilde_eta-R'*Kd*S-M*x(7:9);
end
function [V22,MV2] = bs22(t, x, tau2,dir)
m=116;
Iz=13.1;
X_udot=-167.6;
Y_vdot=-477.2;
N_rdot=-15.9;
Xu=26.9;
Yv=35.8;
Nr=3.5;
Du=241.3;
Dv=503.8;
Dr=76.9;
Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;

% the parameters for the auxiliary control law
M=diag([Mx;My;Mpsi]);
Lambda=diag([1;1;1]);
PI=[0,-1,0;1,0,0;0,0,0];
Kp=diag([1;1;1]);
Kd=diag([1;1;1]);

% the reference trajectory
xr = 0.5*t;
xrdot = 0.5;
xrddot = 0;
xrdddot = 0;

yr = sin(0.5*t+pi/2);
yrdot = 0.5*cos(0.5*t+pi/2);
yrddot = -0.5^2*sin(0.5*t+pi/2);
yrdddot = 0.5^3*cos(0.5*t+pi/2);

thetar = atan2(yrdot,xrdot);
thetardot = (xrdot*yrddot-yrdot*xrddot)/(xrdot^2+yrdot^2);
thetarddot=(1/(xrdot^2+yrdot^2)^2)*((xrddot*yrddot+xrdot*yrdddot-yrddot...
    *xrddot-yrdot*xrdddot)*(xrdot^2+yrdot^2)-(2*xrdot*xrddot+2*yrdot...
    *yrddot)*(xrdot*yrddot-yrdot*xrddot));

eta_r = [xr;yr;thetar]+dir;
eta_r_dot = [xrdot;yrdot;thetardot];
eta_r_ddot = [xrddot;yrddot;thetarddot];

% using backstep technique to contruct the control law
x = x';
eta = x(1:3);
vel = x(4:6);
theta = x(3);
tilde_eta = eta - eta_r;
eta_d_dot = eta_r_dot - Lambda*tilde_eta;

u=vel(1);
v=vel(2);
r=vel(3);
C1=[0;0;My*v];
C2=[0;0;-Mx*u];
C3=[-My*v;Mx*u;0];
C=[C1,C2,C3];
D=diag([-Xu-Du*abs(u);-Yv-Dv*abs(v);-Nr-Dr*abs(r)]);
D=-D;

R = Rot(theta);
eta_dot = R*vel;
v_d = R'*eta_d_dot;
S = eta_dot - eta_d_dot;
Mstar= R*M*R';
Cstar=R*(C-M*R'*vel(3)*R*PI)*R';
Dstar=R*D*R';

thetadot = eta_dot(3);
v_d_dot=-thetadot*R'*PI*eta_d_dot+R'*(eta_r_ddot- Lambda*S);
V22=S'*R*(tau2-M*v_d_dot-C*v_d+D*v_d+R'*Kp*tilde_eta-R'*Kd*S)-S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta;
MV2 = 5*tilde_eta'*Kp*tilde_eta+0.01*S'*Mstar*S;
end

% rotation matrix function
function [ R ] = Rot(theta)
R1=[cos(theta);sin(theta);0];
R2=[-sin(theta);cos(theta);0];
R3=[0;0;1];
R=[R1,R2,R3];
end

