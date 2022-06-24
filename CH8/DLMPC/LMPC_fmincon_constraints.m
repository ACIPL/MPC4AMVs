function [ y,ceq ] = LMPC_fmincon_constraints( u,X0,Hp,Lambda,Kp,Kd,t,T)

ceq=[];

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
% backstepping

M=diag([Mx;My;Mpsi]);
PI=[0,-1,0;1,0,0;0,0,0];
g=0;

psi=X0(3);
r=X0(6);

eta=X0(1:3,1);
vel=X0(4:6,1);
C=CV(vel);
D=DV(vel);

Rpsi=Rot(psi);
Mstar=Rpsi*M*Rpsi';
Cstar=Rpsi*(C-M*Rpsi'*r*Rpsi*PI)*Rpsi';
Dstar=Rpsi*D*Rpsi';
gstar=Rpsi*g;


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

eta_d=[xR;yR;psiR];
eta_d_dot=[xRdot;yRdot;psiRdot];
eta_d_ddot=[xRddot;yRddot;psiRddot];
    
tilde_eta=eta-eta_d;
eta_r_dot=eta_d_dot-Lambda*tilde_eta;

v_r=Rpsi'*eta_r_dot;
    
eta_dot=Rpsi*vel;
S=eta_dot-eta_r_dot;
    
psi_dot=eta_dot(3);
v_r_dot=-psi_dot*Rpsi'*PI*eta_r_dot+Rpsi'*(eta_d_ddot-Lambda*(eta_dot-eta_d_dot));
tau=M*v_r_dot+C*v_r+D*v_r+g-Rpsi'*Kp*tilde_eta-Rpsi'*Kd*S;


V2dot_backstepping=S'*Rpsi*(tau-M*v_r_dot-C*v_r-D*v_r-g+Rpsi'*Kp*tilde_eta)-S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta;

% MPC


nx=length(X0);
nu=length(u); %[U0,U1,..U_N-1]
nu=nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
U=zeros(nu,Hp); 


%==========================================================================
% partition of u | Ui=[Fu_i, Fv_i, Fr_i]'

for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================

%state prediction

Xall=zeros(nx,Hp+1); 
Xplus=X0;
Xall(:,1)=Xplus;

for i=1:1:Hp
    
    Xplus = AUV_Dynamics_discrete(Xplus,U(:,i),T);
    Xall(:,i+1)=Xplus;    
    
end

y=S'*Rpsi*(U(:,1)-M*v_r_dot-C*v_r-D*v_r-g+Rpsi'*Kp*tilde_eta)-S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta-V2dot_backstepping;


end

