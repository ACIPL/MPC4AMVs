function [ y,ceq ] = LMPC_fmincon_constraints( u,X0,Hp,Kp,Kd,t)

ceq=[];

% System coefficients
%==========================================================================
% err_model=0;

B=[0.7974    0.8643    0.8127    0.8270
0.6032    0.5029   -0.5824   -0.5610
0.2945   -0.3302   -0.2847    0.3505];


% m=116;
% Iz=13.1;
% X_udot=-167.6*(1+err_model);
% Y_vdot=-477.2*(1+err_model);
% N_rdot=-15.9*(1+err_model);
% Xu=26.9*(1+err_model);
% Yv=35.8*(1+err_model);
% Nr=3.5*(1+err_model);
% Du=241.3*(1+err_model);
% Dv=503.8*(1+err_model);
% Dr=76.9*(1+err_model);

% Mx=m-X_udot;
% My=m-Y_vdot;
% Mpsi=Iz-N_rdot;
%==========================================================================

% M=diag([Mx;My;Mpsi]);
% PI=[0,-1,0;1,0,0;0,0,0];
% g=0;

psi=X0(3);
% r=X0(6);

eta=X0(1:3,1);
vel=X0(4:6,1);
C=CV(vel);
D=DV(vel);

Rpsi=Rot(psi);

xR=xRef(t);
xRdot=xRefdot(t);   
yR=yRef(t);
yRdot=yRefdot(t);    
psiR=psiRef(xRdot,yRdot);

eta_d=[xR;yR;psiR];
    
tilde_eta=eta-eta_d;

Kd_star=Rpsi'*Kd*Rpsi;
Vdot_PD=-vel'*(D+Kd_star)*vel;

% MPC

% nx=length(X0);
nu=length(u); %[U0,U1,..U_N-1]
nu=nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
U=zeros(nu,Hp); 


%==========================================================================
% partition of u | Ui=[Fu_i, Fv_i, Fr_i]'
%==========================================================================
for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================

Tau=B*U(:,1);

y=vel'*(Tau-C*vel-D*vel+Rpsi'*Kp*tilde_eta)-Vdot_PD;

end

