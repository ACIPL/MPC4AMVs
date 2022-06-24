clc;

eps=1e-1;
stroke=20;% 20 * 10 cm = 2 m
y_scale=0.1; % scale to metric
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

%==========================================================================

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

lb=zeros(nco,1);
ub=10*ones(nco,1);

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

%==========================================================================
for i_sp=1:1:nsp-1
    
    everdone=i_sp*1;
    
    if i_sp == 30
        options = optimset('Algorithm','sqp','Display','iter-detailed','MaxFunEvals',Trail_times);
        ord=5;
        mlti=5;
    end
    
    bbrks=everdone+0:1:everdone+Nbrk;
    X=(everdone+0:1:everdone+Nbrk)';
    Y=zeros(Nbrk+1,1);
    
    for i=1:1:Nbrk+1
        Y(i)=data1_together(1,100*(everdone+i-1)+1);
    end
    %==========================================================================
    knots=augknt(bbrks,mlti);
    nknt=length(knots);
    nco=nknt-ord;
    
    x0=zeros(nco,1);
    x0=x0+1e-5;
    
    lb=zeros(nco,1);
    ub=10*ones(nco,1);
    
    [x,fval] = fmincon(@(x) mycstfun(x,bbrks,mlti,T),x0,[],[],[],[],[],[],@(x) mycons_continuity(x,bbrks,mlti,T,X,Y,stroke,y0,y0d,y02d,eps),options);
    
    ObjVal(i_sp+1,1)=fval;
    
    x=x';
    sp=spmak(knots,x);
    
    
    %==========================================================================
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


%==========================================================================
% Total Plot
T_plt=0.1;
Brk_length=1;
N_every_spline=Brk_length/T_plt;
N_total=nsp*N_every_spline;
Spline_map=zeros(N_total,2);

for i_spl=1:1:nsp
    
    for i_single=1:1:N_every_spline
        
        point_num=(i_spl-1)*N_every_spline+i_single;
        
        Spline_map(point_num,1)=(point_num-1)*T_plt;
        Spline_map(point_num,2)=fnval(SpAry{i_spl,1},Spline_map(point_num,1));
        
    end
    
end

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


figure(3)
plot(Spline_map(:,1),y_scale*Spline_map(:,2));
hold on;

plot(x1(:,1),y1(:,1),'r');
hold on;
plot(x1(:,1),y2(:,1),'r');
grid on;


