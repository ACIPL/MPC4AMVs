function [ fx ] = yaw_cost_2( u,X0,N,Q,Qf,R,P,X_assumed,T)
% u=[Fr0,Fr1,...Fr_N-1]

%==========================================================================
%                    partition weighting matrices
%==========================================================================
q11=Q(1,1);
q22=Q(2,2);
q33=Q(3,3);
q44=Q(4,4);
q55=Q(5,5);
q66=Q(6,6);

Q_omega = diag([q11,q22,q33,q66]);

r11=R(1,1);
r22=R(2,2);
r33=R(3,3);

qf11=Qf(1,1);
qf22=Qf(2,2);
qf33=Qf(3,3);
qf44=Qf(4,4);
qf55=Qf(5,5);
qf66=Qf(6,6);

Qf_omega=diag([qf11,qf22,qf33,qf66]);


%==========================================================================
%               Define the size of needed matrices
%==========================================================================

nu=length(u);
nu=nu/N;

x0=X0(1);
y0=X0(2);
psi0=X0(3);
u0=X0(4);
v0=X0(5);
r0=X0(6);

omega0=[x0;y0;psi0;r0];

n_omega=length(omega0);

Hu = N;
U = zeros(nu,N); 

Omega = zeros(n_omega,N+1);

SurgeVel = zeros(1,N+1);
V = zeros(1,N+1);

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
    V(i) = X_assumed(5,i);
    
end

%==========================================================================
%                         State Prediction
%==========================================================================

Omega(:,1) = omega0;
omega_plus = omega0;

for i=1:1:N
    
    surge_vel = SurgeVel(i);
    v = V(i);
  
    omega_plus =  yaw_subsys_2( omega_plus,surge_vel,v,U(:,i),T );
    Omega(:,i+1) = omega_plus;

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
      
    Omega_r=[xr(i,1);yr(i,1);psi_r(i,1);rr(i,1)];
    
    
    fx = fx + (Omega_r-Omega(:,i))'* Q_omega *(Omega_r-Omega(:,i));

end


for i=1:1:Hu
    
    Ur = Frr(i,1);
    fx = fx + (Ur-U(:,i))' * r33 * (Ur-U(:,i));

end

Omega_r = [xr(N,1);yr(N,1);psi_r(N,1);rr(N,1)];

fx=fx+(Omega_r-Omega(:,N))' * Qf_omega * (Omega_r-Omega(:,N));


end

