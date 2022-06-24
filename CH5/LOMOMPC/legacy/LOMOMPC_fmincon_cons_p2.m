function [ y,ceq ] = LOMOMPC_fmincon_cons_p2( u,J1_star,X0,zeta0,Hp,Q,Qf,R,T,epsi )

nx=length(X0);
nu=length(u); 
nu=nu/Hp; 
U=zeros(nu,Hp); 

%==========================================================================
% partition of u 

for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================
%state prediction

X=zeros(nx,Hp);
Xplus=X0;

ZETA=zeros(2,Hp); % zeta=[s,s_dot];
zeta_plus=zeta0;
y1=zeros(Hp,1);

for i=1:1:Hp
    
    Xplus =AUV_Dynamics_discrete( Xplus,U(1:3,i),T );     
    X(:,i)= Xplus;
    zeta_plus = path_evol_2o( zeta_plus,U(4,i),T );
    ZETA(:,i)=zeta_plus;    
    y1(i)=0-ZETA(2,i);
    
    
end

    xrN = ZETA(1,Hp);
    yrN = sin(ZETA(1,Hp));
    psi_rN = atan(cos(ZETA(1,Hp)));
    
    urN = (px_dot(ZETA(1,Hp))^2+py_dot(ZETA(1,Hp))^2)^(1/2)*ZETA(2,Hp);
    vrN = 0;
    rrN = ZETA(2,Hp)*(px_dot( ZETA(1,Hp) )*py_ddot( ZETA(1,Hp) )-py_dot( ZETA(1,Hp) )*px_ddot( ZETA(1,Hp) ))/(px_dot( ZETA(1,Hp) )^2+py_dot( ZETA(1,Hp) )^2);

    XrN = [xrN,yrN,psi_rN,urN,vrN,rrN]';

    ceq = X(:,Hp)-XrN;


J1 = LOMOMPC_fmincon_cost_p1(u,X0,zeta0,Hp,Q,Qf,R,T);

y2 = J1-J1_star-epsi;

y = [y1;y2];

end

