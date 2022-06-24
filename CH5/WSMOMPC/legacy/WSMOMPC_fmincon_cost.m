function [ fx ] = WSMOMPC_fmincon_cost(u,alpha,uR,X0,s0,Hp,Q1,Qf1,R1,Q2,Qf2,R2,T)

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

nu=length(u); %[U0,U1,..U_N-1]
nu=nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i, v_s]'
nx=length(X0); % X0=[x0,y0,psi0,u0,v0,r0]'

Hu=Hp;
U=zeros(nu,Hp); 
X=zeros(nx,Hp);
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
% partition of u | Ui=[Fu_i, Fv_i, Fr_i, vs_i]'

for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================
% state prediction

Xplus = AUV_Dynamics_discrete( X0,U(1:3,1),T );
X(:,1)=Xplus;
s_plus = path_evol( s0,U(4,1),T );
S(1)=s_plus;


for i=2:1:Hp

    
    Xplus =AUV_Dynamics_discrete( Xplus,U(1:3,i),T );     
    X(:,i)= Xplus;
    
    s_plus = path_evol( s_plus,U(4,i),T );
    S(i)=s_plus;
    

end

%==========================================================================
% reference state 1

for i=1:1:Hp

    xr(i,1)=S(i);
    yr(i,1)=sin(S(i));
    psi_r(i,1)=atan(cos(S(i)));
   
    ur(i,1)=(px_dot(S(i))^2+py_dot(S(i))^2)^(1/2)*U(4,i);
    vr(i,1)=0;
    rr(i,1)=U(4,i)*(px_dot( S(i) )*py_ddot( S(i) )-py_dot( S(i) )*px_ddot( S(i) ))/(px_dot( S(i) )^2+py_dot( S(i) )^2);
        
    Fur(i,1)=-My*vr(i,1)*rr(i,1)+Xu*ur(i,1)+Du*ur(i,1)*abs(ur(i,1));
    Fvr(i,1)= Mx*ur(i,1)*rr(i,1)+Yv*vr(i,1)+Dv*vr(i,1)*abs(vr(i,1));
    Frr(i,1)=(My-Mx)*ur(i,1)*vr(i,1)+Nr*rr(i,1)+Dr*rr(i,1)*abs(rr(i,1))+ ur(i,1)*(-cos(S(i))*(1+cos(S(i))^2)^(3/2)+3*(1+cos(S(i))^2)^(1/2)*sin(S(i))^2*cos(S(i)))/3*(1+cos(S(i))^2)^3;
    vsr(i,1)=ur(i,1)/(1+cos(S(i))^2)^(1/2);
    
    
end    

% reference state 2

[ Xr_N,Sr_N,Ur_N ] = Ref_State_Mat( s0,uR,Hp,T );

%==========================================================================
% cost function

fx=0;

for i=1:1:Hp
      
     Xaug_r1=[xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1);S(i)];
     Xaug=[X(:,i);S(i)]; 
     fx=fx+alpha*(Xaug_r1-Xaug)'*Q1*(Xaug_r1-Xaug);
     
     Xaug_r2=[Xr_N(:,i);Sr_N(i)];
     fx=fx+(1-alpha)*(Xaug_r2-Xaug)'*Q2*(Xaug_r2-Xaug);
     
end

for i=1:1:Hu
    Ur1=[Fur(i,1);Fvr(i,1);Frr(i,1);vsr(i,1)];
    fx=fx+alpha*(Ur1-U(:,i))'*R1*(Ur1-U(:,i));
    
    fx=fx+(1-alpha)*(Ur_N(:,i)-U(:,i))'*R2*(Ur_N(:,i)-U(:,i));
    
end

Xaugf_r1=[xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1);S(Hp)];
Xaugf=[X(:,Hp);S(Hp)]; 
fx=fx+alpha*(Xaugf_r1-Xaugf)'*Qf1*(Xaugf_r1-Xaugf);

Xaugf_r2=[Xr_N(:,Hp);Sr_N(Hp)];
fx=fx+(1-alpha)*(Xaugf_r2-Xaugf)'*Qf2*(Xaugf_r2-Xaugf);


end

