clc;

% System coefficients
%==========================================================================
err_model=0.3;

B=[0.7974    0.8643    0.8127    0.8270
0.6032    0.5029   -0.5824   -0.5610
0.2945   -0.3302   -0.2847    0.3505];

Bplus=pinv(B);
bplus=norm(Bplus,inf);

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
% Backstepping nonlinear control

tau_wave0=[50;50;30];
% tau_wave0=[0;0;0];


M=diag([Mx;My;Mpsi]);
PI=[0,-1,0;1,0,0;0,0,0];
g=0;

Lambda=diag([1;1;1]);
Kp=diag([1;1;1]);
Kd=diag([1;1;1]);

T=0.1;

Tstep=200;

X0=[0.5;0;0;0;0;0];
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
   
%     Backstepping
     tau=M*v_r_dot+C*v_r+D*v_r+g-Rpsi'*Kp*tilde_eta-Rpsi'*Kd*S;
    
    Uall(:,i)=tau;
    
    Thrust_all(:,i)=Bplus*tau;
    tau_wave = tau_wave0;
    tau_wave(1) = tau_wave0(1) * sin(t);
    tau_wave(2) = tau_wave0(2) * cos(0.5*t);
    tau_wave(3) = tau_wave0(3) * sin(0.2*t);
    
    tau=tau+tau_wave;
   
    Xplus  = AUV_Dynamics_discrete_real( Xplus,tau,T );

    Xall(:,i+1)=Xplus;
    
    t=t+T;
end

%==========================================================================
% LMPC

N=5;  % prediction horizon N
Tstep2=Tstep-N; % Total simulation steps in MPC

Xall2=zeros(nx,Tstep2+1);
Uall2=zeros(nu,Tstep2);
Thrust_all2=zeros(n_thrust,Tstep);


Xall2(:,1)=X0;
U0=Uall2(:,1);
Thrust0=Bplus*U0;

%==========================================================================
%                 Limits on Generalized Force and Moments                                 
%==========================================================================
Fu_max=318;
Fv_max=318;
Fr_max=318;
Umax0=[Fu_max;Fv_max;Fr_max];
Umin0=[-Fu_max;-Fv_max;-Fr_max];
Umax=zeros(nu*N,1);
Umin=zeros(nu*N,1);

for i=1:1:N
    
    Umax(1+nu*(i-1):nu*i,1)=Umax0;
    Umin(1+nu*(i-1):nu*i,1)=Umin0;

end

Lbdmax0=[Fu_max;Fv_max];
Lbdmin0=[-Fu_max;-Fv_max];
Lbdmax=zeros(2*N,1);
Lbdmin=zeros(2*N,1);

for i=1:1:N
    
    Lbdmax(1+2*(i-1):2*i,1)=Lbdmax0;
    Lbdmin(1+2*(i-1):2*i,1)=Lbdmin0;

end

Fumax=zeros(N,1);
Fumin=zeros(N,1);

for i=1:1:N
    
    Fumax(1+1*(i-1):1*i,1) = Fu_max;
    Fumin(1+1*(i-1):1*i,1) = -Fu_max;

end

Fvmax=zeros(N,1);
Fvmin=zeros(N,1);

for i=1:1:N
    
    Fvmax(1+1*(i-1):1*i,1) = Fv_max;
    Fvmin(1+1*(i-1):1*i,1) = -Fv_max;

end


Frmax=zeros(N,1);
Frmin=zeros(N,1);

for i=1:1:N
    
    Frmax(1+1*(i-1):1*i,1) = Fr_max;
    Frmin(1+1*(i-1):1*i,1) = -Fr_max;

end

%==========================================================================
%                         Thrust Limits                                  
%==========================================================================
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
%==========================================================================
%                          Weighting Matrices 
%==========================================================================
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
R=diag([r11,r22,r33]);

 R1=0;

qf11=1e3;
qf22=1e3;
qf33=1e2;
qf44=1e2;
qf55=1e2;
qf66=1e2;
Qf=diag([qf11,qf22,qf33,qf44,qf55,qf66]);

%==========================================================================
%                         fmincon setup : DMPC
%==========================================================================
alpha=0.5; 
options = optimset('Algorithm','sqp');

u0_fmincon=Uall(:,1);
for i=2:1:N
    u0_fmincon=[u0_fmincon;Uall(:,i)];
end

Xplus=X0;
P=zeros(nx+nu,N);
U_prev=zeros(nu,1);
R11=diag([0,0,0]);
TIme=zeros(Tstep2,1);
t=0;

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%        Do centralized optimization once for initialization
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    t1=t;
    for j=1:1:N % Generate reference upto N steps ahead
        xR=xRef(t1);
        xRdot=xRefdot(t1);
        xRddot=xRefddot(t1); 
        xRdddot=xRefdddot(t1); 
    
        yR=yRef(t1);
        yRdot=yRefdot(t1);
        yRddot=yRefddot(t1);
        yRdddot=yRefdddot(t1); 
    
        psiR=psiRef(xRdot,yRdot);
        psiRdot=psiRefdot(xRdot,xRddot,yRdot,yRddot);
        psiRddot=psiRefddot(xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot);
        
        uR=sqrt(xRdot^2+yRdot^2);
        vR=0;
        rR=psiRdot;
        
        uRdot=(xRdot^2+yRdot^2)^(-1/2)*(xRdot*xRddot+yRdot*yRddot);
        rRdot=psiRddot;
        
        FuR=Mx*uRdot+Xu*uR+Du*uR*abs(uR);
        FvR=Mx*uR*rR;
        FrR=Mpsi*rRdot+Nr*rR+Dr*rR*abs(rR);
        
        P(:,j)=[xR;yR;psiR;uR;vR;rR;FuR;FvR;FrR];

        t1=t1+T;
    end
    
    rho=(X0-P(1:6,1))'*(X0-P(1:6,1));
    rho=0.1*(rho/2)^(3/2);   
    R12=R1*rho;
    
    u = fmincon(@(u) LMPC_fmincon_cost( u,X0,U_prev,N,Q,Qf,R,R11,P,T),u0_fmincon,[],[],[],[],Umin,Umax,@(u) LMPC_fmincon_constraints( u,Xplus,N,Lambda,Kp,Kd,t,T),options);

    Uini = zeros(nu,N); 
    for i=1:1:N
        for j=1:1:nu
            Uini(j,i)=u((i-1)*nu+j,1);
        end
    end
    Xini = zeros(nx,N+1);
    Xini(:,1) = X0; 
    for i=1:1:N     
            Xplus=AUV_Dynamics_discrete(Xplus,Uini(:,i),T);
            Xini(:,i+1)=Xplus;
    end  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                     Iteration with fmincon
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

X_assumed = Xini;

t3=t;
Xplus3=X0;
tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
u0_surge=tau3(1);
     
for j=2:1:N       
    t3=t3+T;
    Xplus3=AUV_Dynamics_discrete( Xplus3,tau3,T );
    tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
    u0_surge=[u0_surge;tau3(1)];       
end

Xplus = X0;

for i=1:1:Tstep2
    
    tic;      
    t1=t;
    TIme(i,1)=t;
    
    for j=1:1:N % Generate reference upto N steps ahead
        xR=xRef(t1);
        xRdot=xRefdot(t1);
        xRddot=xRefddot(t1); 
        xRdddot=xRefdddot(t1); 
    
        yR=yRef(t1);
        yRdot=yRefdot(t1);
        yRddot=yRefddot(t1);
        yRdddot=yRefdddot(t1); 
    
        psiR=psiRef(xRdot,yRdot);
        psiRdot=psiRefdot(xRdot,xRddot,yRdot,yRddot);
        psiRddot=psiRefddot(xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot);
        
        uR=sqrt(xRdot^2+yRdot^2);
        vR=0;
        rR=psiRdot;
        
        uRdot=(xRdot^2+yRdot^2)^(-1/2)*(xRdot*xRddot+yRdot*yRddot);
        rRdot=psiRddot;
        
        FuR=Mx*uRdot+Xu*uR+Du*uR*abs(uR);
        FvR=Mx*uR*rR;
        FrR=Mpsi*rRdot+Nr*rR+Dr*rR*abs(rR);
        
        P(:,j)=[xR;yR;psiR;uR;vR;rR;FuR;FvR;FrR];

        t1=t1+T;
    end   
    
    u = fmincon(@(u) surge_cost_2( u,Xplus,N,Q,Qf,R,P,X_assumed,T),u0_surge,[],[],[],[],Fumin,Fumax,@(u) surge_constraints( u,Xplus,N,Lambda,Kp,Kd,t,T),options); 
    Fu_pred = u;
    
 %=========================================================================
 %                 Construction of u0_sway
 %=========================================================================

    t3=t;
    Xplus3=Xplus;
    tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
    u0_sway=tau3(2);
    
    for j=2:1:N
        t3=t3+T;
        tau_layer2 = [Fu_pred(j-1);tau3(2);tau3(3)];
        Xplus3=AUV_Dynamics_discrete( Xplus3,tau_layer2,T );
        tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
        u0_sway=[u0_sway;tau3(2)];
    end
     
    u = fmincon(@(u) sway_cost_2( u,Xplus,N,Q,Qf,R,P,X_assumed,T),u0_sway,[],[],[],[],Fvmin,Fvmax,@(u) sway_constraints( u,Xplus,N,Lambda,Kp,Kd,t,Fu_pred,T),options); 
    Fv_pred = u;
 %=========================================================================
 %                 Construction of u0_yaw
 %=========================================================================

    t3=t;
    Xplus3=Xplus;
    tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
    u0_yaw=tau3(3);
    
    for j=2:1:N
        t3=t3+T;
        tau_layer3 = [Fu_pred(j-1);Fv_pred(j-1);tau3(3)];
        Xplus3=AUV_Dynamics_discrete( Xplus3,tau_layer3,T );
        tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
        u0_yaw=[u0_yaw;tau3(3)];
    end
     
    u = fmincon(@(u) yaw_cost_2( u,Xplus,N,Q,Qf,R,P,X_assumed,T),u0_yaw,[],[],[],[],Frmin,Frmax,@(u) yaw_constraints( u,Xplus,N,Lambda,Kp,Kd,t,Fu_pred,Fv_pred,T),options); 
    Fr_pred = u;
    
    tau_actual = [Fu_pred(1);Fv_pred(1);Fr_pred(1)];
    
    Uall2(:,i)=tau_actual;
    Thrust_all2(:,i)=Bplus*tau_actual;
        
    U_prev=tau_actual;
    tau_wave = tau_wave0;
    tau_wave(1) = tau_wave0(1) * sin(t);
    tau_wave(2) = tau_wave0(2) * cos(0.5*t);
    tau_wave(3) = tau_wave0(3) * sin(0.2*t);
    
    tau_actual=tau_actual+tau_wave;
    
    Xplus = AUV_Dynamics_discrete_real( Xplus,tau_actual,T );
    Xall2(:,i+1)=Xplus;    
 %=========================================================================
 %                 Construction of u0_surge 
 %=========================================================================
    t3=t+T;
    Xplus3=Xplus;
    tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
    u0_surge=tau3(1);
     
    for j=2:1:N       
        t3=t3+T;
        Xplus3=AUV_Dynamics_discrete( Xplus3,tau3,T );
        tau3=BS_controller( Xplus3,t3,Lambda,Kp,Kd );
        u0_surge=[u0_surge;tau3(1)];        
    end
 %=========================================================================
 %                    Construction of X_assumed
 %=========================================================================
    X0=Xplus;
    Fu_pred_bar = [Fu_pred(2:N);Fu_pred(N)];
    Fv_pred_bar = [Fv_pred(2:N);Fv_pred(N)];
    Fr_pred_bar = [Fr_pred(2:N);Fr_pred(N)];
    
    Uini = zeros(nu,N); 
    for j=1:1:N     
        Uini(1,j) = Fu_pred_bar(j);
        Uini(2,j) = Fv_pred_bar(j);
        Uini(3,j) = Fr_pred_bar(j);
    end
    Xini = zeros(nx,N+1);
    Xini(:,1) = X0;
    Xplus0 = X0;
    for j=1:1:N     
            Xplus0=AUV_Dynamics_discrete(Xplus0,Uini(:,j),T);
            Xini(:,j+1)=Xplus0;
    end  
    X_assumed = Xini;   
 %=========================================================================
    elapsed_time = toc;

    R11=R12;
    t=t+T;
    i
end

 
%==========================================================================
%                          Plot AUV Trjectory 
%==========================================================================

figure(1)
plot(eta_Ref_all(1,:),eta_Ref_all(2,:),'k','LineWidth',2);
hold on;
plot(Xall(1,1:Tstep2),Xall(2,1:Tstep2),'b','LineWidth',2);
hold on;
plot(Xall2(1,1:Tstep2),Xall2(2,1:Tstep2),'r','LineWidth',2);
hold on;

%==========================================================================
%                        Plot Generalized Force
%==========================================================================
figure(2)
subplot(3,1,1)
% plot(TIme(1:Tstep2,1),Uall(1,1:Tstep2),'b');
% hold on;
plot(TIme(1:Tstep2,1),Uall2(1,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(3,1,2)
% plot(TIme(1:Tstep2,1),Uall(2,1:Tstep2),'b');
% hold on;
plot(TIme(1:Tstep2,1),Uall2(2,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(3,1,3)
% plot(TIme(1:Tstep2,1),Uall(3,1:Tstep2),'b');
% hold on;
plot(TIme(1:Tstep2,1),Uall2(3,1:Tstep2),'r','LineWidth',2);
hold on;
%==========================================================================
%                        Plot Thrusts
%==========================================================================
figure(3)
subplot(4,1,1)
% plot(TIme(1:Tstep2,1),Thrust_all(1,1:Tstep2),'b');
% hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(1,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(4,1,2)
% plot(TIme(1:Tstep2,1),Thrust_all(2,1:Tstep2),'b');
% hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(2,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(4,1,3)
% plot(TIme(1:Tstep2,1),Thrust_all(3,1:Tstep2),'b');
% hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(3,1:Tstep2),'r','LineWidth',2);
hold on;
 
subplot(4,1,4)
% plot(TIme(1:Tstep2,1),Thrust_all(4,1:Tstep2),'b');
% hold on;
plot(TIme(1:Tstep2,1),Thrust_all2(4,1:Tstep2),'r','LineWidth',2);
hold on;
 
%==========================================================================
%                     Plot Reference vs. AUV x,y psi
%==========================================================================
figure(4)
subplot(3,1,1)
% plot(TIme(1:Tstep2,1),eta_Ref_all(1,1:Tstep2),'k');
% hold on;
% plot(TIme(1:Tstep2,1),Xall(1,1:Tstep2),'b');
% hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(1,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(3,1,2)
% plot(TIme(1:Tstep2,1),eta_Ref_all(2,1:Tstep2),'k');
% hold on;
% plot(TIme(1:Tstep2,1),Xall(2,1:Tstep2),'b');
% hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(2,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(3,1,3)
% plot(TIme(1:Tstep2,1),eta_Ref_all(3,1:Tstep2),'k');
% hold on;
% plot(TIme(1:Tstep2,1),Xall(3,1:Tstep2),'b');
% hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(3,1:Tstep2),'r','LineWidth',2);
hold on;

%==========================================================================
%                        Plot AUV u,v,r
%==========================================================================
figure(5)
subplot(3,1,1)
% plot(TIme(1:Tstep2,1),Xall(4,1:Tstep2),'b');
% hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(4,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(3,1,2)
% plot(TIme(1:Tstep2,1),Xall(5,1:Tstep2),'b');
% hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(5,1:Tstep2),'r','LineWidth',2);
hold on;

subplot(3,1,3)
% plot(TIme(1:Tstep2,1),Xall(6,1:Tstep2),'b');
% hold on;
grid on;
plot(TIme(1:Tstep2,1),Xall2(6,1:Tstep2),'r','LineWidth',2);
hold on;


%==========================================================================
%                    Calculate MSE on Each Reference
%==========================================================================
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

    MSEx1=MSEx1/Tstep2
    MSEx2=MSEx2/Tstep2
    
    MSEy1=MSEy1/Tstep2
    MSEy2=MSEy2/Tstep2
    
    MSEpsi1=MSEpsi1/Tstep2
    MSEpsi2=MSEpsi2/Tstep2


