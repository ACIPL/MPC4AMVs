function [ cost ] = LOMOMPC_fmincon_cost_p1(u,X0,zeta0,Hp,Q,Qf,R,T)

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
q11 = Q(1,1);
q22 = Q(2,2);
q33 = Q(3,3);
q44 = Q(4,4);
q55 = Q(5,5);
q66 = Q(6,6);
q77 = Q(7,7);

r11 = R(1,1);
r22 = R(2,2);
r33 = R(3,3);
r44 = R(4,4);

qf11 = Qf(1,1);
qf22 = Qf(2,2);
qf33 = Qf(3,3);
qf44 = Qf(4,4);
qf55 = Qf(5,5);
qf66 = Qf(6,6);
qf77 = Qf(7,7);
%==========================================================================

nu=length(u); 
nu=nu/Hp;
nx=length(X0); 

Hu=Hp;
U=zeros(nu,Hp); 
X=zeros(nx,Hp);
ZETA=zeros(2,Hp); % zeta=[s,s_dot];


xr=zeros(Hp,1);
yr=zeros(Hp,1);
psi_r=zeros(Hp,1);
ur=zeros(Hp,1);
vr=zeros(Hp,1);
rr=zeros(Hp,1);

P=zeros(9,Hp); % P(:,i) = [xri;yri;psi_ri;uri;vri;rri;urdi;vrdi;rrdi]; 

zeta_r=zeros(Hp,2);

ur_d=zeros(Hp,1);
rr_d=zeros(Hp,1);

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


for i=2:1:Hp

    
    Xplus =AUV_Dynamics_discrete( Xplus,U(1:3,i),T );     
    X(:,i)= Xplus;
    
    zeta_plus = path_evol_2o( zeta_plus,U(4,i),T );
    ZETA(:,i)=zeta_plus;


end

%==========================================================================
% reference state

for i=1:1:Hp
    
    xr(i,1)=ZETA(1,i);
    yr(i,1)=sin(ZETA(1,i));
    psi_r(i,1)=atan(cos(ZETA(1,i)));
    
    ur(i,1)=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)^(1/2)*ZETA(2,i);
    vr(i,1)=0;
    rr(i,1)=ZETA(2,i)*(px_dot( ZETA(1,i) )*py_ddot( ZETA(1,i) )-py_dot( ZETA(1,i) )*px_ddot( ZETA(1,i) ))/(px_dot( ZETA(1,i) )^2+py_dot( ZETA(1,i) )^2);
    
    vr_d = 0;
   
    urd1=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)^(-1/2)*(px_dot(ZETA(1,i))*px_ddot(ZETA(1,i))+py_dot(ZETA(1,i))*py_ddot(ZETA(1,i)))*U(4,i);
    urd2=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)^(1/2)*U(4,i);
    ur_d(i,1)=urd1+urd2;

    rrd1=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)*(px_dot(ZETA(1,i))*py_dddot(ZETA(1,i))-py_dot(ZETA(1,i))*px_dddot(ZETA(1,i)))-2*(px_dot(ZETA(1,i))*py_ddot(ZETA(1,i))-py_dot(ZETA(1,i))*px_ddot(ZETA(1,i)))*(px_dot(ZETA(1,i))*px_ddot(ZETA(1,i))+py_dot(ZETA(1,i))*py_ddot(ZETA(1,i)));
    rrd2=px_dot(ZETA(1,i))*py_ddot(ZETA(1,i))-py_dot(ZETA(1,i))*px_ddot(ZETA(1,i));
    deno=px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2;
    rr_d(i,1)=(rrd1/deno^2+rrd2/deno)*U(4,i);
    
    P(:,i) = [xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1);ur_d(i,1);vr_d;rr_d(i,1)];
    
    Fur(i,1)=Mx*ur_d(i,1)+Xu*ur(i,1)+Du*ur(i,1)*abs(ur(i,1));
    Fvr(i,1)= Mx*ur(i,1)*rr(i,1);
    Frr(i,1)=Mpsi*rr_d(i,1)+Nr*rr(i,1)+Dr*rr(i,1)*abs(rr(i,1));
    vsr(i,1)=U(4,i); 
   
end    

%==========================================================================
% cost function

 cost=0;
 Q_bar = diag([q11,q22,q33,q44,q55,q66]);
 R_bar = diag([r11,r22,r33]);
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
     cost = cost+Xe'*Q_bar*Xe+Te'*R_bar*Te;
     
%      fx = fx+q77*(ZETA(1,i)-ZETA(1,i))^2 + r44*(U(4,i)-vsr(i,1))^2; % equal to zero
end


 Qf_bar = diag([qf11,qf22,qf33,qf44,qf55,qf66]);

 Xr=[xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
 Xe = ErrorState( X(:,Hp), Xr );
 cost = cost+Xe'*Qf_bar*Xe;
 
%  fx = fx+qf77*(ZETA(1,Hp)-ZETA(1,Hp))^2; % equal to zero



end

