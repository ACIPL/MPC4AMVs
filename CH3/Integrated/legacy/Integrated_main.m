clear; 
clc;

load('SensorData.mat')

eps=1e-1;
stroke=20; % 20 * 10 cm = 2 m
y_scale=0.1;
T=0.1;
nsp=40;
Nbrk=3;
ord=4;
mlti=4;

SpAry=cell(nsp,1);
ObjVal=zeros(nsp,1);
ObjVal_real=zeros(nsp,1);
Cost_tatal=0;
Trail_times=1200;

everdone=0;
bbrks=everdone+0:1:everdone+Nbrk;
X=(everdone+0:1:everdone+Nbrk)';
Y=zeros(Nbrk+1,1);

for i=1:1:Nbrk+1
Y(i)=data1_together(1,100*(everdone+i-1)+1);
end

knots=augknt(bbrks,mlti);
nknt=length(knots);
nco=nknt-ord;

x0=zeros(nco,1);
x0=x0+1e-5;

options = optimset('Algorithm','interior-point','Display','iter-detailed','MaxFunEvals',Trail_times);

[x,fval] = fmincon(@(x) mycstfun(x,bbrks,mlti,T),x0,[],[],[],[],[],[],@(x) mycons(x,bbrks,mlti,T,X,Y,stroke),options);

ObjVal(1,1)=fval;
x=x';
sp=spmak(knots,x);

y0= fnval(sp,everdone+1);
y0d= fnval(fnder(sp,1),everdone+1);
y02d= fnval(fnder(sp,2),everdone+1);
SpAry(1,1)={sp}; 

fprime=fnder(sp);
f2prime=fnder(sp,2);
f3prime=fnder(sp,3);

fx=0;

for x=bbrks(1):T:bbrks(2)

kprimex=(fnval(f3prime,x)*(1+fnval(fprime,x)^2)^(3/2)-fnval(f2prime,x)*(3/2)*(1+fnval(fprime,x)^2)^(1/2)*2*fnval(fprime,x)*fnval(f2prime,x))/(1+fnval(fprime,x)^2)^3;
kprimex=kprimex^2;
f0x=kprimex/(1+fnval(fprime,x)^2)^(1/2);
fx=fx+T*f0x;

end

ObjVal_real(1,1)=fx;

Cost_tatal=sum(ObjVal_real);
iter_round=1; 
everdone=1;

M=10;
T_plt=T;
Brk_length=1;
N_every_spline=Brk_length/T_plt;
N_total=nsp*N_every_spline;
Spline_map2=zeros(N_total,5);% rescale the data into metric (in meter) | [xr;S(xr);S'(xr);S"(xr);S"'(xr)];
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

N = 8;  % prediction horizon N
Nx = m-N;
Tstep = Nx-N-1;
Tstep2 = Tstep-N; 

sp_required = ceil(N/M);

X0=[0;1.5;0;0;0;0];
nx=length(X0);
nu=3;

Xall=zeros(nx,Tstep+1);
Uall=zeros(nu,Tstep);
eta_Ref_all=zeros(3,Tstep);

Xall(:,1)=X0;

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

u0_fmincon=Uall(:,1);
for i=2:1:N
    u0_fmincon=[u0_fmincon;Uall(:,i)];
end

Xplus=X0;
P=zeros(nx+nu,N);

options = optimset('Algorithm','interior-point','Display','iter-detailed','MaxFunEvals',Trail_times);
if sp_required >1
    for i_sp=1:1:(sp_required-1)
        everdone=i_sp*1;
        
        bbrks=everdone+0:1:everdone+Nbrk;
        X=(everdone+0:1:everdone+Nbrk)';
        Y=zeros(Nbrk+1,1);
        
        for i=1:1:Nbrk+1
            Y(i)=data1_together(1,100*(everdone+i-1)+1);
        end
        
        knots=augknt(bbrks,mlti);
        nknt=length(knots);
        nco=nknt-ord;
        
        x0=zeros(nco,1);
        x0=x0+1e-5;
        
        [x,fval] = fmincon(@(x) mycstfun(x,bbrks,mlti,T),x0,[],[],[],[],[],[],@(x) mycons_continuity(x,bbrks,mlti,T,X,Y,stroke,y0,y0d,y02d,eps),options);
        ObjVal(i_sp+1,1)=fval;
        
        x=x';
        sp=spmak(knots,x);
        
        y0= fnval(sp,everdone+1);
        y0d= fnval(fnder(sp,1),everdone+1);
        y02d= fnval(fnder(sp,2),everdone+1);

        SpAry(i_sp+1,1)={sp};
        
        fprime=fnder(sp);
        f2prime=fnder(sp,2);
        f3prime=fnder(sp,3);

        fx=0;
        
        for x=bbrks(1):T:bbrks(2)
            kprimex=(fnval(f3prime,x)*(1+fnval(fprime,x)^2)^(3/2)-fnval(f2prime,x)*(3/2)*(1+fnval(fprime,x)^2)^(1/2)*2*fnval(fprime,x)*fnval(f2prime,x))/(1+fnval(fprime,x)^2)^3;
            kprimex=kprimex^2;
            f0x=kprimex/(1+fnval(fprime,x)^2)^(1/2);
            fx=fx+T*f0x;
        end
        ord=4;
        mlti=4; 
        ObjVal_real(i_sp,1)=fx;
        Cost_tatal=sum(ObjVal_real);
        iter_round=iter_round+1;
    end
end

%==========================================================================
t=0; Tstep2=80;
for i=1:1:Tstep2

    if mod(i+N-1,M) == 0
        
        options = optimset('Algorithm','interior-point','Display','iter-detailed','MaxFunEvals',Trail_times);
        everdone=iter_round*1;
        if everdone == 30
            options = optimset('Algorithm','sqp','Display','iter-detailed','MaxFunEvals',Trail_times);
            ord=5;
            mlti=5;
        elseif everdone >= 31
            options = optimset('Algorithm','sqp','Display','iter-detailed','MaxFunEvals',Trail_times);
        end
        
        bbrks=everdone+0:1:everdone+Nbrk;
        X=(everdone+0:1:everdone+Nbrk)';
        Y=zeros(Nbrk+1,1);

        for l=1:1:Nbrk+1
            Y(l)=data1_together(1,100*(everdone+l-1)+1);
        end
        
        knots=augknt(bbrks,mlti);
        nknt=length(knots);
        nco=nknt-ord;

        x0=zeros(nco,1);
        x0=x0+1e-5;
        
        [x,fval] = fmincon(@(x) mycstfun(x,bbrks,mlti,T),x0,[],[],[],[],[],[],@(x) mycons_continuity(x,bbrks,mlti,T,X,Y,stroke,y0,y0d,y02d,eps),options);

        ObjVal(everdone+1,1)=fval;

        x=x';
        sp=spmak(knots,x);
        
        y0= fnval(sp,everdone+1);
        y0d= fnval(fnder(sp,1),everdone+1);
        y02d= fnval(fnder(sp,2),everdone+1);
        SpAry(everdone+1,1)={sp};
        
        fprime=fnder(sp);
        f2prime=fnder(sp,2);
        f3prime=fnder(sp,3);
        
        fx=0;
        for x=bbrks(1):T:bbrks(2)
            kprimex=(fnval(f3prime,x)*(1+fnval(fprime,x)^2)^(3/2)-fnval(f2prime,x)*(3/2)*(1+fnval(fprime,x)^2)^(1/2)*2*fnval(fprime,x)*fnval(f2prime,x))/(1+fnval(fprime,x)^2)^3;
            kprimex=kprimex^2;
            f0x=kprimex/(1+fnval(fprime,x)^2)^(1/2);
            fx=fx+T*f0x;
        end
        ord=4;
        mlti=4; 
        
        ObjVal_real(everdone,1)=fx;
        Cost_tatal=sum(ObjVal_real);
        
        iter_round=iter_round+1;
    end
    
    xR=xRef(t); 
    xRdot=xRefdot(t);
    
    i_spl=floor((i-1)/M)+1;
    Spline_map2(i,1)=(i-1)*T_plt;
    Spline_map2(i,2)=y_scale*fnval(SpAry{i_spl,1},Spline_map2(i,1));
    Spline_map2(i,3)=y_scale*fnval(fnder(SpAry{i_spl,1},1),Spline_map2(i,1)); 
    Spline_map2(i,4)=y_scale*fnval(fnder(SpAry{i_spl,1},2),Spline_map2(i,1));
    Spline_map2(i,5)=y_scale*fnval(fnder(SpAry{i_spl,1},3),Spline_map2(i,1)); 
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
       
        i_spl=floor((i+j-1)/M)+1;
        Spline_map2(i+j,1)=(i+j-1)*T_plt;
        Spline_map2(i+j,2)=y_scale*fnval(SpAry{i_spl,1},Spline_map2(i+j,1));
        Spline_map2(i+j,3)=y_scale*fnval(fnder(SpAry{i_spl,1},1),Spline_map2(i+j,1)); 
        Spline_map2(i+j,4)=y_scale*fnval(fnder(SpAry{i_spl,1},2),Spline_map2(i+j,1));
        Spline_map2(i+j,5)=y_scale*fnval(fnder(SpAry{i_spl,1},3),Spline_map2(i+j,1)); 
       
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
    options = optimset('Algorithm','sqp','Display','iter-detailed');
    
    tic;
    u = fmincon(@(u) NMPC_fmincon_cost( u,X0,N,Q,Qf,R,P,T),u0_fmincon,[],[],[],[],-Umax,Umax,@(u) NMPC_fmincon_constraints( u,Xplus,N,P,T),options);
    toc;
    
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

