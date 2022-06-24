clc;

% System coefficients
%==========================================================================
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
%==========================================================================

tau_wave=[0;0;0];


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

Xall=zeros(nx,Tstep+1);
Uall=zeros(nu,Tstep);
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
    
    tau=M*v_r_dot+C*v_r+D*v_r+g-Rpsi'*Kp*tilde_eta-Rpsi'*Kd*S;
   
    Uall(:,i)=tau;
    tau=tau+tau_wave;
    
    Xplus  = AUV_Dynamics_discrete_real( Xplus,tau,T );

    Xall(:,i+1)=Xplus;
    
    t=t+T;
end

figure(1)
plot(eta_Ref_all(1,:),eta_Ref_all(2,:),'r');
hold on;
plot(Xall(1,1:Tstep),Xall(2,1:Tstep));
hold on;

figure(2)
subplot(3,1,1)
plot(Uall(1,:),'b');

subplot(3,1,2)
plot(Uall(2,:),'b');

subplot(3,1,3)
plot(Uall(3,:),'b');

hold on;

figure(3)
subplot(3,1,1)
plot(Xall(1,:),'b');
grid on;

subplot(3,1,2)
plot(Xall(2,:),'b');
grid on;

subplot(3,1,3)
plot(Xall(3,:),'b');
grid on;

figure(4)
subplot(3,1,1)
plot(Xall(4,:),'b');
grid on;

subplot(3,1,2)
plot(Xall(5,:),'b');
grid on;

subplot(3,1,3)
plot(Xall(6,:),'b');
grid on;

