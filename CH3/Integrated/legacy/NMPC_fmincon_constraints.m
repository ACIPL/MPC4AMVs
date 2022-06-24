function [ y,ceq ] = NMPC_fmincon_constraints( u,X0,Hp,P,T)

ceq=[];

%==========================================================================
% MPC

nx=length(X0);
nu=length(u); %[U0,U1,..U_N-1]
nu=nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
U=zeros(nu,Hp); 

xr=zeros(Hp,1);
yr=zeros(Hp,1);
psi_r=zeros(Hp,1);
ur=zeros(Hp,1);
vr=zeros(Hp,1);
rr=zeros(Hp,1);
urdot=zeros(Hp,1);
vrdot=zeros(Hp,1);
rrdot=zeros(Hp,1);


%==========================================================================
% partition of u | Ui=[Fu_i, Fv_i, Fr_i]'

for i=1:1:Hp
    for j=1:1:nu

        U(j,i)=u((i-1)*nu+j,1);
         
    end
end

%==========================================================================

% reference state

for i=1:1:Hp

    xr(i,1)=P(1,i);
    yr(i,1)=P(2,i);
    psi_r(i,1)=P(3,i);
    ur(i,1)=P(4,i);
    vr(i,1)=P(5,i);
    rr(i,1)=P(6,i);
    
    urdot(i,1)=P(7,i);
    vrdot(i,1)=P(8,i);
    rrdot(i,1)=P(9,i);
    
end    

%==========================================================================


%state prediction

Xall=zeros(nx,Hp+1); 
Xplus=X0;
Xall(:,1)=Xplus;

for i=1:1:Hp-1
    
    Xplus = AUV_Dynamics_discrete(Xplus,U(:,i),T);
    Xall(:,i+1)=Xplus;    
    
end

XN = Xall(:,Hp);
XrN = [ xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
XeN = ErrorState( XN, XrN );
xeN = XeN(1);
yeN = XeN(2);
psi_eN = XeN(3);
ueN = XeN(4);
veN = XeN(5);
reN = XeN(6);

y1 = abs(xeN) - abs(ueN);
y2 = abs(yeN) - abs(veN);
y3 = abs(psi_eN) - abs(reN);
y4 = -xeN*ueN;
y5 = -yeN*veN;
y6 = -psi_eN*reN;

y=[y1;y2;y3;y4;y5;y6];


end

