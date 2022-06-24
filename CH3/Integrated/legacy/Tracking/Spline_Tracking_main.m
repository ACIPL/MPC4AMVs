clc;

%==========================================================================
% Prepare reference data
T=0.1;
y_scale=0.1;

nsp=40;
stroke=20; % 20 * 10 cm = 2 m
T_plt=T;

Brk_length=1;
N_every_spline=Brk_length/T_plt;
N_total=nsp*N_every_spline;
Spline_map2=zeros(N_total,5);% rescale the data into metric (in meter) | [xr;S(xr);S'(xr);S"(xr);S"'(xr)];

for i_spl=1:1:nsp
    for i_single=1:1:N_every_spline
        
        point_num=(i_spl-1)*N_every_spline+i_single;
        
        Spline_map2(point_num,1)=(point_num-1)*T_plt;
        Spline_map2(point_num,2)=y_scale*fnval(SpAry{i_spl,1},Spline_map2(point_num,1));
        Spline_map2(point_num,3)=y_scale*fnval(fnder(SpAry{i_spl,1},1),Spline_map2(point_num,1)); 
        Spline_map2(point_num,4)=y_scale*fnval(fnder(SpAry{i_spl,1},2),Spline_map2(point_num,1));
        Spline_map2(point_num,5)=y_scale*fnval(fnder(SpAry{i_spl,1},3),Spline_map2(point_num,1));
        
    end
    
end
[m,n]=size(Spline_map2);

X=(0:1:nsp)';
Y=X;

for i=1:1:nsp+1
Y(i)=data1_together(1,100*(i-1)+1);
end

x1=X(1):T_plt:X(nsp+1);
x1=x1';
nx1=length(x1);

y1=zeros(nx1,1);
y2=y1;
for i=1:1:nx1

y1(i)=y_scale*plygfun1(x1(i),X,Y);
y2(i)=y_scale*plygfun2(x1(i),X,Y,stroke);

end

%==========================================================================
N=8;  % prediction horizon N
Nx = m-N;
Tstep = Nx-N-1;
Tstep2=Tstep-N; 

X0=[0;1.5;0;0;0;0];
nx=length(X0);
nu=3;

Xall=zeros(nx,Tstep+1);
Uall=zeros(nu,Tstep);
eta_Ref_all=zeros(3,Tstep);

Xall(:,1)=X0;
%==========================================================================
% NMPC error dynamics

Xall2=zeros(nx,Tstep2+1);
Uall2=zeros(nu,Tstep2);

Xall2(:,1)=X0;
U0=Uall2(:,1:N);

Fu_max=2000;
Fv_max=2000;
Fr_max=1000;

Umax0=[Fu_max;Fv_max;Fr_max];
Umin0=[-Fu_max;-Fv_max;-Fr_max];
Umax=zeros(nu*N,1);
Umin=zeros(nu*N,1);

for i=1:1:N
    
    Umax(1+nu*(i-1):nu*i,1)=Umax0;
    Umin(1+nu*(i-1):nu*i,1)=Umin0;

end

% Weighting Matrices 
Q = 0.4 * eye(6);

R = 0.01 * eye(3);

Qf = 0.5 * eye(6);

K1 = 0.5 * eye(3);
K2 = -1 * eye(3);
K = [K1 K2];

options = optimset('Algorithm','sqp');

u0_fmincon=Uall(:,1);
for i=2:1:N
    u0_fmincon=[u0_fmincon;Uall(:,i)];
end

Xplus=X0;
t=0;

P=zeros(nx+nu,N);

tic;
for i=1:1:Tstep2
    
    xR=xRef(t); 
    xRdot=xRefdot(t);
    yR=Spline_map2(i,2);
    yRdot=Spline_map2(i,3);
    
    psiR=psiRef(xRdot,yRdot);
    eta_Ref_all(:,i)=[xR;yR;psiR];
    t1=t;
    for j=1:1:N % Generate reference upto N steps ahead
        xR=xRef(t1);
        xRdot=xRefdot(t1);
        xRddot=xRefddot(t1); 
        xRdddot=xRefdddot(t1); 
       
        yR=Spline_map2(i+j-1,2);
        yRdot=Spline_map2(i+j-1,3);
        yRddot=Spline_map2(i+j-1,4);
        yRdddot=Spline_map2(i+j-1,5);

        psiR=psiRef(xRdot,yRdot);
        psiRdot=psiRefdot(xRdot,xRddot,yRdot,yRddot);
        psiRddot=psiRefddot(xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot);
 
        uR=sqrt(xRdot^2+yRdot^2);
        vR=0;
        rR=psiRdot;
        
        uRdot=(xRdot^2+yRdot^2)^(-1/2)*(xRdot*xRddot+yRdot*yRddot);
        vRdot=0;
        rRdot=psiRddot;

        P(:,j)=[xR;yR;psiR;uR;vR;rR;uRdot;vRdot;rRdot];

        t1=t1+T;
    end
    
    u = fmincon(@(u) NMPC_fmincon_cost( u,X0,N,Q,Qf,R,P,T),u0_fmincon,[],[],[],[],-Umax,Umax,@(u) NMPC_fmincon_constraints( u,Xplus,N,P,T),options);
    
    u_actual=u(1:nu,1);
    Uall2(:,i)=u_actual;
    
 %==========================================================================
 % Construction of a feasible control u0_fmincon using terminal controller 
 %==========================================================================
    t2=t;
    Xplus2=Xplus;
   
    for k=1:1:N  % partition of u
        for j=1:1:nu

            U(j,k)=u((k-1)*nu+j,1);

        end
    end
    
    for j=1:1:N-1 % state prediction
        
        t2=t2+T;
        Xplus2=AUV_Dynamics_discrete( Xplus2,U(:,j),T );
        
    end
    
    XN = Xplus2;
    XrN = P(1:6,N);
    XeN = ErrorState( XN, XrN );
    u0_fmincon=[u(nu+1:nu*N,1);K*XeN];
 %==========================================================================      
     
    Xplus = AUV_Dynamics_discrete( Xplus,Uall2(:,i),T );
    Xall2(:,i+1)=Xplus;

    X0=Xplus;
    
    i
    
    t=t+T;
end
   toc;  



%==========================================================================
% Plot
figure(1)
plot(x1(:,1),y1(:,1),'g','LineWidth',2);
hold on;
plot(x1(:,1),y2(:,1),'g','LineWidth',2);
grid on;
plot(eta_Ref_all(1,1:Tstep2),eta_Ref_all(2,1:Tstep2),'k','LineWidth',2);
hold on;
plot(Xall2(1,1:Tstep2),Xall2(2,1:Tstep2),'r','LineWidth',2);
hold on;

figure(2)
subplot(3,1,1)
plot(Uall2(1,:),'r');
hold on;

subplot(3,1,2)
plot(Uall2(2,:),'r');
hold on;

subplot(3,1,3)
plot(Uall2(3,:),'r');
hold on;

