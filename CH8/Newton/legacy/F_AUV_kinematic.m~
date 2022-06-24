function [ F ] = F_AUV_kinematic( X0,U,P,N,Q,R,Qf,T,Umax,tau )
% U has size of [9,N]; U=[U0,U1,...,U_N-1]
% Ui=[Fu;Fv;Fr;u4;u5;u6;mu1;mu2;mu3]

% Umax=[Fu_max;Fv_max;Fr_max]

% P has size of [6,N+1]; P=[Xd0,Xd1,...,Xd_N]
% Xdi=[xd;yd;psi_d;ud;vd;rd]

% X evolution throuth system dynamics
%==========================================================================
[nU,mU]=size(U);

m=length(X0);
X=zeros(m,N+1); % X=[X0,X1,...,XN]

X(:,1)=X0;
for j=1:1:N
   X(:,j+1)=( X(:,j),U(:,j),T );
end
%==========================================================================
%Lambda propagation backwardly
%==========================================================================
Lambda=zeros(m,N+1); % Lambda=[lambda0,lambda1,...,lambdaN]
Lambda(:,N+1)=PHIx_AUV_kinematic( X(:,N+1),Qf,P(:,N+1));


for i=N:-1:1
    Hx = Hx_AUV_kinematic( X(:,i),Lambda(:,i+1),U,Q,P(:,i));
    Lambda(:,i)=Lambda(:,i+1)+Hx*T;
end
%==========================================================================

F=zeros(nU*N,1);

for k=1:1:N
    F1=Hu_AUV_kinematic( X(:,k),Lambda(:,k+1),U(:,k),R,tau,Umax);
    F(nU*(k-1)+1:nU*k,1)=F1;
end
end

