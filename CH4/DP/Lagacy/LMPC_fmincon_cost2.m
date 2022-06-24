function [ fx ] = LMPC_fmincon_cost2( u,X0,Hp,Q,Qf,R,P,T)

nu=length(u);
nu=nu/Hp;
nx=length(X0);

Hu=Hp;
U=zeros(nu,Hp); 
X=zeros(nx,Hp);

xr=zeros(Hp,1);
yr=zeros(Hp,1);
psi_r=zeros(Hp,1);
ur=zeros(Hp,1);
vr=zeros(Hp,1);
rr=zeros(Hp,1);
Fur=zeros(Hp,1);
Fvr=zeros(Hp,1);
Frr=zeros(Hp,1);

T1=zeros(Hp,1);
T2=zeros(Hp,1);
T3=zeros(Hp,1);
T4=zeros(Hp,1);

B=[0.7974    0.8643    0.8127    0.8270
0.6032    0.5029   -0.5824   -0.5610
0.2945   -0.3302   -0.2847    0.3505];

%==========================================================================
% partition of u

for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================
% state prediction

Xplus = AUV_Dynamics_discrete( X0,B*U(:,1),T );

X(:,1)=Xplus;

for i=2:1:Hp

    
    Xplus =AUV_Dynamics_discrete( Xplus,B*U(:,i),T );
        
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
    
    Fur(i,1)=P(7,i);
    Fvr(i,1)=P(8,i);
    Frr(i,1)=P(9,i);
    
    Tau=[Fur(i,1);Fvr(i,1);Frr(i,1)];
    Thrust=pinv(B)*Tau;
    
    T1(i,1)=Thrust(1);
    T2(i,1)=Thrust(2);
    T3(i,1)=Thrust(3);
    T4(i,1)=Thrust(4);
    
end    

%==========================================================================
% cost function

fx=0;

for i=1:1:Hp
    
    Xr=[xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1)];
   
    fx=fx+(Xr-X(:,i))'*Q*(Xr-X(:,i));

end

for i=1:1:Hu
    
    fx=fx+U(:,i)'*R*U(:,i);
end

Xr=[xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
fx=fx+(Xr-X(:,Hp))'*Qf*(Xr-X(:,Hp));

end
