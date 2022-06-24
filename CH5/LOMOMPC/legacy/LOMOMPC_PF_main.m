clc;
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
M=100; % total simulation steps

N=10;  
T=0.1;
nx=6; 
nu=4; 
uR=1;

epsi=15;
epsi_max=15;
bias=5;

Fu_max=500;
Fv_max=500;
Fr_max=500;
v_smax=100;
v_smin=-100;
Umax0=[Fu_max;Fv_max;Fr_max;v_smax];
Umin0=[-Fu_max;-Fv_max;-Fr_max;v_smin];
Umax=zeros(nu*N,1);
Umin=zeros(nu*N,1);

for i=1:1:N
    
    Umax(1+nu*(i-1):nu*i,1)=Umax0;
    Umin(1+nu*(i-1):nu*i,1)=Umin0;

end

%==========================================================================
% Weighting Matrices for Priority 1
q11=1e5;
q22=1e5;
q33=1e3;
q44=1e3;
q55=1e-3;
q66=1e-3;
q77=1e-3;
Q=diag([q11,q22,q33,q44,q55,q66,q77]);

r11=1e-3;
r22=1e-3;
r33=1e-3;
r44=1e-3;
R=diag([r11,r22,r33,r44]);

qf11=100;
qf22=100;
qf33=10;
qf44=100;
qf55=1e-3;
qf66=1e-3;
qf77=1e-3;
Qf=diag([qf11,qf22,qf33,qf44,qf55,qf66,qf77]);

%==========================================================================
% Weighting Matrices for Priority 2
q11_p2=1e4;
q22_p2=1e4;
q33_p2=1e2;
q44_p2=1e3;
q55_p2=1e-3;
q66_p2=1e-3;
q77_p2=1e3;
Q_prio2=diag([q11_p2,q22_p2,q33_p2,q44_p2,q55_p2,q66_p2,q77_p2]);


r11_p2=1e-3;
r22_p2=1e-3;
r33_p2=1e-3;
r44_p2=1e-3;
R_prio2=diag([r11_p2,r22_p2,r33_p2,r44_p2]);

qf11_p2=100;
qf22_p2=100;
qf33_p2=10;
qf44_p2=100;
qf55_p2=1e-3;
qf66_p2=1e-3;
qf77_p2=1e2;
Qf_prio2=diag([qf11_p2,qf22_p2,qf33_p2,qf44_p2,qf55_p2,qf66_p2,qf77_p2]);

Qepsi = 1e-1*eye(7,7);
%==========================================================================

X=zeros(nx,N+1); % X=[X0,X1,...,XN]
U=zeros(nu,N);   % U=[U0,U1,...,U_N-1]

Xall=zeros(nx,M+1);
Uall=zeros(nu,M);
ZETAall=zeros(2,M+1);
Sall=zeros(1,M+1);

J1all=zeros(1,M);
J2all=zeros(1,M);

X0=[0.5;0;0;0;0;0];
zeta0 = [0;0];
s0 = zeta0(1,1);

Xall(:,1) = X0;
ZETAall(:,1) = zeta0;
Sall(1) = s0;

U0=Uall(:,1:N); 

options = optimset('Algorithm','sqp');

u0_fmincon=zeros(nu*N,1);
Xplus = X0;
zeta_plus = zeta0;

for i=1:1:M
    
    tic  
    [u,J1_star] = fmincon(@(u) LOMOMPC_fmincon_cost_p1( u,Xplus,zeta_plus,N,Q,Qf,R,T),u0_fmincon,[],[],[],[],Umin,Umax,@(u) LOMOMPC_fmincon_cons_p1( u,Xplus,zeta_plus,N,T),options);
    toc
    J1all(i)=J1_star;
    u1=u;
    
    u0_fmincon_p2=u1;
    
    tic  
    [u,J2_star] = fmincon(@(u) LOMOMPC_fmincon_cost_p2( u,uR,Xplus,zeta_plus,N,Q_prio2,Qf_prio2,R_prio2,T),u0_fmincon_p2,[],[],[],[],Umin,Umax,@(u) LOMOMPC_fmincon_cons_p2( u,J1_star,Xplus,zeta_plus,N,Q,Qf,R,T,epsi ),options);
    toc
    J2all(i) = J2_star;
   
    u_actual = u(1:nu,1);
    u_actual_AUV = u(1:nu-1,1);
    Uall(:,i) = u_actual;
    u0_fmincon = Uini( u,Xplus,zeta_plus,N,T );

    Xplus = AUV_Dynamics_discrete( Xplus,u_actual_AUV,T );
    Xall(:,i+1)=Xplus;
    zeta_plus = path_evol_2o( zeta_plus,u(4,1),T );
    ZETAall(:,i+1)=zeta_plus;

    s_plus = zeta_plus(1,1);
    Sall(i+1) = s_plus;
       
    xr = ZETAall(1,i+1);
    yr = sin(ZETAall(1,i+1));
    psi_r = atan(cos(ZETAall(1,i+1)));
    ur = (px_dot(ZETAall(1,i+1))^2+py_dot(ZETAall(1,i+1))^2)^(1/2)*ZETAall(2,i+1);
    vr = 0;
    rr = ZETAall(2,i+1)*(px_dot( ZETAall(1,i+1) )*py_ddot( ZETAall(1,i+1) )-py_dot( ZETAall(1,i+1) )*px_ddot( ZETAall(1,i+1) ))/(px_dot( ZETAall(1,i+1) )^2+py_dot( ZETAall(1,i+1) )^2);
    Xr = [xr;yr;psi_r;ur;vr;rr];
    Xe = ErrorState( Xall(:,i+1), Xr );    
    Xaug_tilde = [Xe;0];
    
    xr_prev = ZETAall(1,i);
    yr_prev = sin(ZETAall(1,i));
    psi_r_prev = atan(cos(ZETAall(1,i)));
    ur_prev = (px_dot(ZETAall(1,i))^2+py_dot(ZETAall(1,i))^2)^(1/2)*ZETAall(2,i);
    vr_prev = 0;
    rr_prev = ZETAall(2,i)*(px_dot( ZETAall(1,i) )*py_ddot( ZETAall(1,i) )-py_dot( ZETAall(1,i) )*px_ddot( ZETAall(1,i) ))/(px_dot( ZETAall(1,i) )^2+py_dot( ZETAall(1,i) )^2);
    vr_d_prev = 0;
    urd1_prev = (px_dot(ZETAall(1,i))^2+py_dot(ZETAall(1,i))^2)^(-1/2)*(px_dot(ZETAall(1,i))*px_ddot(ZETAall(1,i))+py_dot(ZETAall(1,i))*py_ddot(ZETAall(1,i)))*Uall(4,i);
    urd2_prev = (px_dot(ZETAall(1,i))^2+py_dot(ZETAall(1,i))^2)^(1/2)*Uall(4,i);
    ur_d_prev = urd1_prev + urd2_prev;
    rrd1_prev = (px_dot(ZETAall(1,i))^2+py_dot(ZETAall(1,i))^2)*(px_dot(ZETAall(1,i))*py_dddot(ZETAall(1,i))-py_dot(ZETAall(1,i))*px_dddot(ZETAall(1,i)))-2*(px_dot(ZETAall(1,i))*py_ddot(ZETAall(1,i))-py_dot(ZETAall(1,i))*px_ddot(ZETAall(1,i)))*(px_dot(ZETAall(1,i))*px_ddot(ZETAall(1,i))+py_dot(ZETAall(1,i))*py_ddot(ZETAall(1,i)));
    rrd2_prev = px_dot(ZETAall(1,i))*py_ddot(ZETAall(1,i))-py_dot(ZETAall(1,i))*px_ddot(ZETAall(1,i));
    deno_prev = px_dot(ZETAall(1,i))^2+py_dot(ZETAall(1,i))^2;
    rr_d_prev = (rrd1_prev/deno_prev^2+rrd2_prev/deno_prev)*Uall(4,i);
    P_prev = [xr_prev;yr_prev;psi_r_prev;ur_prev;vr_prev;rr_prev;ur_d_prev;vr_d_prev;rr_d_prev];
    Xr_prev = [xr_prev;yr_prev;psi_r_prev;ur_prev;vr_prev;rr_prev];
    Xe_prev = ErrorState( Xall(:,i), Xr_prev );
    f1_prev = F1( Xall(:,i),Xe_prev,P_prev );
    f2_prev = F2( Xall(:,i),Xe_prev,P_prev );
    f3_prev = F3( Xall(:,i),Xe_prev,P_prev );
    tau_u_prev = (Uall(1,i)-f1_prev)/Mx;  
    tau_v_prev = (Uall(2,i)-f2_prev)/My;
    tau_r_prev = (Uall(3,i)-f3_prev)/Mpsi;
    Te_prev = [tau_u_prev;tau_v_prev;tau_r_prev];

    Xaug_prev_tilde = [Xe_prev;0];
    Uaug_prev_tilde = [Te_prev;0];
    
    epsi0 = Xaug_prev_tilde'*Q*Xaug_prev_tilde + Uaug_prev_tilde'*R*Uaug_prev_tilde - Xaug_tilde'*Qepsi*Xaug_tilde;    
    if epsi0 > epsi_max
        epsi1 = epsi_max;
    else
        epsi1 = epsi0;
    end    
    epsi = epsi1 + bias;
        
    i
    
end
%==========================================================================
xs=Sall;
ys=Sall;
for i=1:1:M+1
    
    ys(i)=sin(xs(i));
    
end

OT=zeros(M,1);

for i=1:1:M
   
    OT(i)=i*T;
    
end

%plot figures

figure(1)
plot(Xall(1,:),Xall(2,:),'k:');
hold on;


figure(2)
subplot(2,2,1)
plot(OT(:,1),Uall(1,:),'b--'), title('Fu')
hold on;
grid on;

subplot(2,2,2)
plot(OT(:,1),Uall(2,:),'b--'), title('Fv')
hold on;
grid on;

subplot(2,2,3)
plot(OT(:,1),Uall(3,:),'b--'), title('Fr')
hold on;
grid on;

subplot(2,2,4)
plot(OT(:,1),Uall(4,:),'b--'), title('vs')
hold on;
grid on;


Uall(:,1:10)

figure(3)
plot(Xall(1,:),Xall(2,:),'b');
hold on;
grid on;

figure(4)
plot(OT(:,1),Xall(4,1:M),'k:');
title('Surge Velocity of the Vehicle')
ylabel('u [m/s]')
xlabel('Time [sec.]')
hold on;
grid on;

figure(5)

subplot(2,1,1)
plot(OT(:,1),J1all(1:M),'b--');
title('Value Function of J1')
ylabel('J1^*')
xlabel('Time [sec.]')
hold on;
grid on;

subplot(2,1,2)
plot(OT(:,1),J2all(1:M),'b--');
title('Value Function of J2')
ylabel('J2^*')
xlabel('Time [sec.]')
hold on;
grid on;


