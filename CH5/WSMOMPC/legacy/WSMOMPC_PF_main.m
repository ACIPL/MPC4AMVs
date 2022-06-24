clc;
%==========================================================================
M=100; % total simulation steps

N=10;  % prediction horizon N
T=0.1; % sampling period T
nx=6; % X=[x;y;psi;u;v;r]
nu=4; % Ui=[Fu;Fv;Fr,vs]
uR=1;


Fu_max=500;
Fv_max=500;
Fr_max=500;
v_smax=100;
v_smin=0;
Umax0=[Fu_max;Fv_max;Fr_max;v_smax];
Umin0=[-Fu_max;-Fv_max;-Fr_max;v_smin];

Umax = zeros(nu*N,1);
Umin = zeros(nu*N,1);

for i=1:1:N
    
    Umax(1+nu*(i-1):nu*i,1)=Umax0;
    Umin(1+nu*(i-1):nu*i,1)=Umin0;

end


alpha=0.99;
beta=1e0;

Q1 = diag([1e5,1e5,1e2,1e-1,1e-1,1e-1,1e-1]);
R1 = diag([1e-3,1e-3,1e-3,1e-3]);
Qf1 = diag([1e2,1e2,1e1,1e-3,1e-3,1e-3,1e-3]);

Q2 = diag([1e0,1e0,1e0,1e3,1e-1,1e-1,1e-1]);
R2 = diag([1e-3,1e-3,1e-3,1e-3]);
Qf2 = diag([1e-1,1e-1,1e-1,1e2,1e-3,1e-3,1e-3]);

%==========================================================================

X=zeros(nx,N+1); % X=[X0,X1,...,XN]
U=zeros(nu,N);   % U=[U0,U1,...,U_N-1] | Ui=[Fu_i, Fv_i, Fr_i, s_i]'

Xall=zeros(nx,M+1);
Uall=zeros(nu,M);
Sall=zeros(1,M+1);

X0=[0.5;0;0;0;0;0];
s0=0;
Xall(:,1)=X0;
Sall(1)=s0;

U0=Uall(:,1:N); % U0=0

options = optimset('Algorithm','sqp');

u0_fmincon=zeros(nu*N,1);
Xplus=X0;
s_plus=s0;
xR=s0;
yR=sin(s0);
psi_R=atan(cos(s0));
for i=1:1:M
     
    err=(Xall(1,i)-xR)^2+(Xall(2,i)-yR)^2+(Xall(3,i)-psi_R)^2;
    alpha=alpha_logistic( err,beta );
           
    tic  
    u = fmincon(@(u) WSMOMPC_fmincon_cost( u,alpha,uR,Xplus,s_plus,N,Q1,Qf1,R1,Q2,Qf2,R2,T),u0_fmincon,[],[],[],[],Umin,Umax,[],options);
    toc
    
    u_actual=u(1:nu,1);
    u_actual_AUV=u(1:nu-1,1);
    Uall(:,i)=u_actual;
    
    Xplus = AUV_Dynamics_discrete( Xplus,u_actual_AUV,T );
    Xall(:,i+1)=Xplus;
    s_plus = path_evol( s_plus,u(4,1),T );
    Sall(i+1)=s_plus;

    u0_fmincon=[u(nu+1:nu*N,1);u(nu*(N-1)+1:nu*N,1)];
    
    xR=s_plus;
    yR=sin(s_plus);
    psi_R=atan(cos(s_plus));
    
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

% plot figures
figure(1)
xd=0:0.1:10; % size of 101
yd=sin(xd);  % size of 101
Pall=[xd;yd]';
plot(Pall(:,1),Pall(:,2),'g')
hold on;
plot(xs,ys,'r')
hold on;
plot(Xall(1,:),Xall(2,:),'b');
hold on;

figure(2)
subplot(2,2,1)
plot(OT(:,1),Uall(1,:),'b'), title('Fu')
hold on;
grid on;

subplot(2,2,2)
plot(OT(:,1),Uall(2,:),'b'), title('Fv')
hold on;
grid on;

subplot(2,2,3)
plot(OT(:,1),Uall(3,:),'b'), title('Fr')
hold on;
grid on;

subplot(2,2,4)
plot(OT(:,1),Uall(4,:),'b'), title('vs')
hold on;
grid on;


Uall(:,1:10)

figure(3)
plot(Xall(1,:),Xall(2,:),'b');
hold on;
grid on;

figure(4)
plot(OT(:,1),Xall(4,1:M));
title('Surge Velocity of the Vehicle')
ylabel('u')
hold on;
grid on;



