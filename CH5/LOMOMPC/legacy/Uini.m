function [ u_initial ] = Uini( u,X0,zeta0,N,T )
%zeta=[s,s_dot]'

nu=length(u); 
nu=nu/N; 
nx=length(X0);


U=zeros(nu,N); 
X=zeros(nx,N);
ZETA=zeros(2,N);

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

%==========================================================================
% partition of u 

for i=1:1:N
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================
% state prediction

Xplus = AUV_Dynamics_discrete( X0,U(1:3,1),T );
X(:,1)=Xplus;
zeta_plus = path_evol( zeta0,U(4,1),T );
ZETA(:,1)=zeta_plus;


for i=2:1:N

    
    Xplus =AUV_Dynamics_discrete( Xplus,U(1:3,i),T );     
    X(:,i)= Xplus;
    
    zeta_plus = path_evol( zeta_plus,U(4,i),T );
    ZETA(:,i)=zeta_plus;
   
end

%==========================================================================
% local controller
    ur=(px_dot(ZETA(1,N))^2+py_dot(ZETA(1,N))^2)^(1/2)*ZETA(2,N);
    vr=0;
    rr=ZETA(2,N)*(px_dot( ZETA(1,N) )*py_ddot( ZETA(1,N) )-py_dot( ZETA(1,N) )*px_ddot( ZETA(1,N) ))/(px_dot( ZETA(1,N) )^2+py_dot( ZETA(1,N) )^2);
    ur_d=0;
    rr_d=0;
    
    Fur=Mx*ur_d+Xu*ur+Du*ur*abs(ur);
    Fvr= Mx*ur*rr;
    Frr=Mpsi*rr_d+Nr*rr+Dr*rr*abs(rr);

 u_initial=[u(nu+1:nu*N,1);[Fur;Fvr;Frr;0]];

end

