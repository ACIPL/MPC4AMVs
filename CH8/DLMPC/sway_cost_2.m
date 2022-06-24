function [ fx ] = sway_cost_2( u,X0,N,Q,Qf,R,P,X_assumed,T)

%==========================================================================
%                    Partition weighting matrices
%==========================================================================

q11=Q(1,1);
q22=Q(2,2);
q33=Q(3,3);
q44=Q(4,4);
q55=Q(5,5);
q66=Q(6,6);

Q_zeta = diag([q11,q22,q33,q55]);

r11=R(1,1);
r22=R(2,2);
r33=R(3,3);

qf11=Qf(1,1);
qf22=Qf(2,2);
qf33=Qf(3,3);
qf44=Qf(4,4);
qf55=Qf(5,5);
qf66=Qf(6,6);

Qf_zeta = diag([qf11,qf22,qf33,qf55]);

%==========================================================================
%                Define the size of needed matrices
%==========================================================================
nu=length(u);
nu=nu/N;

x0=X0(1);
y0=X0(2);
psi0=X0(3);
u0=X0(4);
v0=X0(5);
r0=X0(6);

zeta0=[x0;y0;psi0;v0];

n_zeta=length(zeta0);

Hu=N;
U=zeros(nu,N); 

Zeta = zeros(n_zeta,N+1);

SurgeVel = zeros(1,N+1);
YawRate = zeros(1,N+1);

xr=zeros(N+1,1);
yr=zeros(N+1,1);
psi_r=zeros(N+1,1);
ur=zeros(N+1,1);
vr=zeros(N+1,1);
rr=zeros(N+1,1);
Fur=zeros(N+1,1);
Fvr=zeros(N+1,1);
Frr=zeros(N+1,1);

%==========================================================================
%                        Partition of u
%==========================================================================

for i=1:1:N
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================
%                         Assumed State 
%==========================================================================

for i=1:1:N

    SurgeVel(i) = X_assumed(4,i);
    YawRate(i) = X_assumed(6,i);
    
end

%==========================================================================
%                         State Prediction
%==========================================================================

Zeta(:,1) = zeta0;
zeta_plus = zeta0;
for i=1:1:N
    
    surge_vel = SurgeVel(i);
    r = YawRate(i);
    zeta_plus =  sway_subsys_2( zeta_plus,surge_vel,r,U(:,i),T );
    Zeta(:,i+1) = zeta_plus;

end

%==========================================================================
%                         Reference State
%==========================================================================

for i=1:1:N
    
    xr(i,1)=P(1,i);
    yr(i,1)=P(2,i);
    psi_r(i,1)=P(3,i);
    ur(i,1)=P(4,i);
    vr(i,1)=P(5,i);
    rr(i,1)=P(6,i);
    Fur(i,1)=P(7,i);
    Fvr(i,1)=P(8,i);
    Frr(i,1)=P(9,i);
    
end    

%==========================================================================
%                          Cost function
%==========================================================================

fx=0;

for i=1:1:N
     
    Zeta_r=[xr(i,1);yr(i,1);psi_r(i,1);vr(i,1)];
    
    fx=fx+(Zeta_r-Zeta(:,i))'*Q_zeta*(Zeta_r-Zeta(:,i));

end


for i=1:1:Hu
    
    Ur = Fvr(i,1);
    fx = fx + (Ur-U(:,i))'* r22 * (Ur-U(:,i));

end

Zeta_r = [xr(N,1);yr(N,1);psi_r(N,1);vr(N,1)];

fx = fx+(Zeta_r-Zeta(:,N))'* Qf_zeta *(Zeta_r-Zeta(:,N));


end

