function [ cost ] = LOMOMPC_fmincon_cost_p2(u,uR,X0,zeta0,Hp,Q_prio2,Qf_prio2,R_prio2,T)

% System coefficients
%==========================================================================
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
%==========================================================================

s0 = zeta0(1,1);

nu=length(u);
nu=nu/Hp; 
nx=length(X0); 
Hu=Hp;
U=zeros(nu,Hp); 
X=zeros(nx,Hp);
ZETA=zeros(2,Hp);
S=zeros(1,Hp);

xr=zeros(Hp,1);
yr=zeros(Hp,1);
psi_r=zeros(Hp,1);
ur=zeros(Hp,1);
vr=zeros(Hp,1);
rr=zeros(Hp,1);
sr=zeros(Hp,1);

Fur=zeros(Hp,1);
Fvr=zeros(Hp,1);
Frr=zeros(Hp,1);
vsr=zeros(Hp,1);

%==========================================================================
% partition of u

for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================
% state prediction

Xplus = AUV_Dynamics_discrete( X0,U(1:3,1),T );
X(:,1)=Xplus;
zeta_plus = path_evol_2o( zeta0,U(4,1),T );
ZETA(:,1)=zeta_plus;
S(1)=ZETA(1,1);

for i=2:1:Hp
   
    Xplus =AUV_Dynamics_discrete( Xplus,U(1:3,i),T );     
    X(:,i)= Xplus;
    
    zeta_plus = path_evol_2o( zeta_plus,U(4,i),T );
    ZETA(:,i) = zeta_plus;
    S(i) = ZETA(1,i);
end

%==========================================================================
% reference state generation

[ Xr_N,Sr_N,Ur_N ] = Ref_State_Mat( s0,uR,Hp,T );


%==========================================================================
% cost function

cost=0;

for i=1:1:Hp

     Xaug=[X(:,i);S(i)];   
     Xaug_r=[Xr_N(:,i);Sr_N(i)];
     cost=cost+(Xaug_r-Xaug)'*Q_prio2*(Xaug_r-Xaug);
     
end

for i=1:1:Hu
    
    cost=cost+(Ur_N(:,i)-U(:,i))'*R_prio2*(Ur_N(:,i)-U(:,i));

end

    Xaugf=[X(:,Hp);S(Hp)];
    Xaugf_r=[Xr_N(:,Hp);Sr_N(Hp)];
    cost=cost+(Xaugf_r-Xaugf)'*Qf_prio2*(Xaugf_r-Xaugf);
end

