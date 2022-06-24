function [ fx ] = NMPC_fmincon_cost( u,X0,Hp,Q,Qf,R,P,T)

% System coefficients
%==========================================================================
err_model=0;

m=116;
Iz=13.1;
X_udot=-100.6*(1+err_model);
Y_vdot=-477.2*(1+err_model);
N_rdot=-15.9*(1+err_model);
Xu=26.9*(1+err_model);
Yv=35.8*(1+err_model);
Nr=3.5*(1+err_model);
Du=100.3*(1+err_model);
Dv=503.8*(1+err_model);
Dr=76.9*(1+err_model);

Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;

%==========================================================================

nu=length(u);
nu=nu/Hp;
nx=length(X0);

U=zeros(nu,Hp); 
X=zeros(nx,Hp);

xr=zeros(Hp,1);
yr=zeros(Hp,1);
psi_r=zeros(Hp,1);
ur=zeros(Hp,1);
vr=zeros(Hp,1);
rr=zeros(Hp,1);
urdot=zeros(Hp,1);
vrdot=zeros(Hp,1);
rrdot=zeros(Hp,1);


%==========================================================================
% partition of u

for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================
% state prediction

Xplus =AUV_Dynamics_discrete( X0,U(:,1),T );

X(:,1)=Xplus;

for i=2:1:Hp

    
    Xplus =AUV_Dynamics_discrete( Xplus,U(:,i),T );
        
    X(:,i)= Xplus;

end

%==========================================================================
% reference state

for i=1:1:Hp

    xr(i,1)=P(1,i);
    yr(i,1)=P(2,i);
    psi_r(i,1)=P(3,i);
    ur(i,1)=P(4,i);
    vr(i,1)=P(5,i);
    rr(i,1)=P(6,i);
    
    urdot(i,1)=P(7,i);
    vrdot(i,1)=P(8,i);
    rrdot(i,1)=P(9,i);
    
end    

%==========================================================================
% cost function

fx=0;

for i=1:1:Hp
    
    Xr=[xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1)];
    Xe = ErrorState( X(:,i), Xr );
    f1 = F1( X(:,i),Xe,P(:,i) );
    f2 = F2( X(:,i),Xe,P(:,i) );
    f3 = F3( X(:,i),Xe,P(:,i) );
    tau_u = (U(1,i)-f1)/Mx;  
    tau_v = (U(2,i)-f2)/My;
    tau_r = (U(3,i)-f3)/Mpsi;
    Te = [tau_u;tau_v;tau_r];
    fx=fx+Xe'*Q*Xe+Te'*R*Te;

end

Xr=[xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
Xe = ErrorState( X(:,Hp), Xr );
fx=fx+Xe'*Qf*Xe;

end
