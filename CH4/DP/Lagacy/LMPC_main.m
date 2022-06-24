clc;

% System coefficients
%==========================================================================
err_model=0;

B=[0.7974    0.8643    0.8127    0.8270
0.6032    0.5029   -0.5824   -0.5610
0.2945   -0.3302   -0.2847    0.3505];

Bplus=pinv(B);

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
%==========================================================================
% Nonlinear PD control

% tau_wave=[4;4;0];
tau_wave=[0;0;0];


M=diag([Mx;My;Mpsi]);
PI=[0,-1,0;1,0,0;0,0,0];
g=0;

Lambda=diag([1;1;1]);
Kp=diag([10;10;10]);
Kd=diag([10;10;10]);

T=0.1;

Tstep=200;

X0=[5;5;-pi/2;0;0;0];

nx=length(X0);
nu=3;
n_thrust=4;

Xall=zeros(nx,Tstep+1);
Uall=zeros(nu,Tstep);
Thrust_all=zeros(n_thrust,Tstep);

eta_Ref_all=zeros(3,Tstep);

t=0;

Xall(:,1)=X0;
Xplus=X0;
for i=1:1:Tstep
 
    xR=xRef(t);
    xRdot=xRefdot(t);
    xRddot=xRefddot(t); 
    xRdddot=xRefdddot(t); 
    
    yR=yRef(t);
    yRdot=yRefdot(t);
    yRddot=yRefddot(t);
    yRdddot=yRefdddot(t); 
    
    psiR=psiRef(xRdot,yRdot);
    psiRdot=psiRefdot(xRdot,xRddot,yRdot,yRddot);
    psiRddot=psiRefddot(xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot);
    
    eta_Ref_all(:,i)=[xR;yR;psiR];
    
    eta=Xall(1:3,i);
    vel=Xall(4:6,i);
    psi=Xall(3,i);
    
    eta_d=[xR;yR;psiR];
    eta_d_dot=[xRdot;yRdot;psiRdot];
    eta_d_ddot=[xRddot;yRddot;psiRddot];
    
    tilde_eta=eta-eta_d;
    eta_r_dot=eta_d_dot-Lambda*tilde_eta;
    
    Rpsi=Rot(psi);
    v_r=Rpsi'*eta_r_dot;
    
    eta_dot=Rpsi*vel;
    S=eta_dot-eta_r_dot;
    
    psi_dot=eta_dot(3);
    v_r_dot=-psi_dot*Rpsi'*PI*eta_r_dot+Rpsi'*(eta_d_ddot-Lambda*(eta_dot-eta_d_dot));
    
    C=CV(vel);
    D=DV(vel);
    
%     PD control
    tau=-Rpsi'*(Kp*tilde_eta+Kd*eta_dot);
    
    Uall(:,i)=tau;
    
    Thrust_all(:,i)=Bplus*tau;
    
    tau=tau+tau_wave;
   
    Xplus  = AUV_Dynamics_discrete_real( Xplus,tau,T );

    Xall(:,i+1)=Xplus;
    
    t=t+T;
end

%==========================================================================
% LMPC

N = 3;  % prediction horizon N
Tstep2=Tstep - N; 

Xall2=zeros(nx,Tstep2+1);
Uall2=zeros(nu,Tstep2);
Thrust_all2=zeros(n_thrust,Tstep);


Xall2(:,1)=X0;
U0=Uall2(:,1);
Thrust0=Bplus*U0;

T1_max=500;
T2_max=500;
T3_max=500;
T4_max=500;
Thrust_max0=[T1_max;T2_max;T3_max;T4_max];
Thrust_min0=[-T1_max;-T2_max;-T3_max;-T4_max];
Thrust_max=zeros(n_thrust*N,1);
Thrust_min=zeros(n_thrust*N,1);

for i=1:1:N
    
    Thrust_max(1+n_thrust*(i-1):n_thrust*i,1)=Thrust_max0;
    Thrust_min(1+n_thrust*(i-1):n_thrust*i,1)=Thrust_min0; 

end

% Weighting Matrices 
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

r1_11=3*1e-2;
r1_22=3*1e-2;
r1_33=3*1e-2;
r1_44=3*1e-2;
R1=diag([r1_11,r1_22,r1_33,r1_44]);

qf11=1e3;
qf22=1e3;
qf33=1e2;
qf44=1e1;
qf55=1e1;
qf66=1e1;
Qf=diag([qf11,qf22,qf33,qf44,qf55,qf66]);

options = optimset('Algorithm','sqp');

u0_fmincon = zeros(n_thrust*N,1);

Xplus=X0;
t=0;

P=zeros(nx+nu,N);

U_prev=zeros(n_thrust,1);

R11=diag([0,0,0,0]);

TIme=zeros(Tstep2,1);

tic;
for i=1:1:Tstep2
        
    t1=t;
    TIme(i,1)=t;
    
    for j=1:1:N % Generate reference upto N steps ahead
        xR=xRef(t1);
        xRdot=xRefdot(t1);
%         xRddot=xRefddot(t1); 
%         xRdddot=xRefdddot(t1); 
    
        yR=yRef(t1);
        yRdot=yRefdot(t1);
%         yRddot=yRefddot(t1);
%         yRdddot=yRefdddot(t1); 
    
        psiR=psiRef(xRdot,yRdot);
        psiRdot=psiRefdot(xRdot,xRddot,yRdot,yRddot);
        psiRddot=psiRefddot(xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot);
        
        uR=sqrt(xRdot^2+yRdot^2);
        vR=0;
        rR=psiRdot;
        
        uRdot=0;
        rRdot=psiRddot;
        
        FuR=Mx*uRdot+Xu*uR+Du*uR*abs(uR);
        FvR=Mx*uR*rR;
        FrR=Mpsi*rRdot+Nr*rR+Dr*rR*abs(rR);
        
        P(:,j)=[xR;yR;psiR;uR;vR;rR;FuR;FvR;FrR];

        t1=t1+T;
    end
    
%     rho=(X0-P(1:6,1))'*(X0-P(1:6,1));
%     rho=0.1*(rho/2)^(3/2);   
%     
%     R12=R1;
    
%     u = fmincon(@(u) LMPC_fmincon_cost( u,X0,U_prev,N,Q,Qf,R,R11,P,T),u0_fmincon,[],[],[],[],Thrust_min,Thrust_max,@(u) LMPC_fmincon_constraints( u,Xplus,N,Lambda,Kp,Kd,t,T),options);
    u = fmincon(@(u) LMPC_fmincon_cost2( u,X0,N,Q,Qf,R,P,T),u0_fmincon,[],[],[],[],Thrust_min,Thrust_max,@(u) LMPC_fmincon_constraints( u,Xplus,N,Kp,Kd,t),options);
    
    u_actual=u(1:n_thrust,1);
    Thrust_all2(:,i)=u_actual;
    
    tau_actual=B*u_actual;
    Uall2(:,i)=tau_actual;
    
%     U_prev=u_actual;
    
    tau_actual=tau_actual+tau_wave;
    
    Xplus = AUV_Dynamics_discrete_real( Xplus,tau_actual,T );
    Xall2(:,i+1)=Xplus;

 %==========================================================================
 % Construction of a feasible control u0_fmincon using PD controller 
 
    t2=t+T;
    Xplus2=Xplus;
    tau2=PD_controller( Xplus2,t2,Lambda,Kp,Kd );
    u0_fmincon=Bplus*tau2;
    
    for j=2:1:N
        
        t2=t2+T;
        Xplus2=AUV_Dynamics_discrete( Xplus2,tau2,T );
        tau2=PD_controller( Xplus2,t2,Lambda,Kp,Kd );
        u0_fmincon=[u0_fmincon;Bplus*tau2];
        
    end
    u0_fmincon
 %==========================================================================   
    
    X0=Xplus;
    
    i
  
%     R11=R12;
    t=t+T;
end
   toc;  

%==========================================================================
% Plot

figure(1)
plot(Xall(1,1:Tstep2),Xall(2,1:Tstep2),'b');
hold on;
plot(Xall2(1,1:Tstep2),Xall2(2,1:Tstep2),'r');
hold on;
grid on;

figure(2)
subplot(3,1,1)
plot(TIme(1:Tstep2,1),Uall(1,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Uall2(1,1:Tstep2),'r');


subplot(3,1,2)
plot(TIme(1:Tstep2,1),Uall(2,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Uall2(2,1:Tstep2),'r');


subplot(3,1,3)
plot(TIme(1:Tstep2,1),Uall(3,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Uall2(3,1:Tstep2),'r');


figure(3)
subplot(4,1,1)
plot(TIme(1:Tstep2,1),Thrust_all(1,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(1,1:Tstep2),'r');


subplot(4,1,2)
plot(TIme(1:Tstep2,1),Thrust_all(2,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(2,1:Tstep2),'r');


subplot(4,1,3)
plot(TIme(1:Tstep2,1),Thrust_all(3,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(3,1:Tstep2),'r');

subplot(4,1,4)
plot(TIme(1:Tstep2,1),Thrust_all(4,1:Tstep2),'b');
hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(4,1:Tstep2),'r');


figure(4)
subplot(3,1,1)
plot(TIme(1:Tstep2,1),Xall(1,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(1,1:Tstep2),'r');


subplot(3,1,2)
plot(TIme(1:Tstep2,1),Xall(2,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(2,1:Tstep2),'r');


subplot(3,1,3)
plot(TIme(1:Tstep2,1),Xall(3,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(3,1:Tstep2),'r');



figure(5)
subplot(3,1,1)
plot(TIme(1:Tstep2,1),Xall(4,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(4,1:Tstep2),'r');


subplot(3,1,2)
plot(TIme(1:Tstep2,1),Xall(5,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(5,1:Tstep2),'r');


subplot(3,1,3)
plot(TIme(1:Tstep2,1),Xall(6,1:Tstep2),'b');
hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(6,1:Tstep2),'r');


