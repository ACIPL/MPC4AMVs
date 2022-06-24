function [ Hx ] = Hx_AUV_kinematic( X,lambda,U,Q,p )

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

x=X(1);
y=X(2);
psi=X(3);
u=X(4);
v=X(5);
r=X(6);

xd=p(1);
yd=p(2);
psi_d=p(3);
ud=p(4);
vd=p(5);
rd=p(6);


lambda1=lambda(1);
lambda2=lambda(2);
lambda3=lambda(3);
lambda4=lambda(4);
lambda5=lambda(5);
lambda6=lambda(6);

Fu=U(1);
Fv=U(2);
Fr=U(3);

q11=Q(1,1);
q22=Q(2,2);
q33=Q(3,3);
q44=Q(4,4);
q55=Q(5,5);
q66=Q(6,6);

Hx1=q11*(x-xd);
Hx2=q22*(y-yd);
Hx3=q33*(psi-psi_d)-lambda1*u*sin(psi)-lambda1*v*cos(psi)+lambda2*u*cos(psi)-lambda2*v*sin(psi);
Hx4=q44*(u-ud)+lambda1*cos(psi)+lambda2*sin(psi)-lambda4*(Xu/Mx)-lambda5*(Mx/My)*r+lambda6*((Mx-My)/Mpsi)*v-2*lambda4*(Du/Mx)*abs(u);
Hx5=q55*(v-vd)-lambda1*sin(psi)+lambda2*cos(psi)+lambda4*(My/Mx)*r-lambda5*(Yv/My)+lambda6*((Mx-My)/Mpsi)*u-2*lambda5*(Dv/My)*abs(v);
Hx6=q66*(r-rd)+lambda3+lambda4*(My/Mx)*v-lambda5*(Mx/My)*u-lambda6*(Nr/Mpsi)-2*lambda6*(Dr/Mpsi)*abs(r);


Hx=[Hx1;Hx2;Hx3;Hx4;Hx5;Hx6];

end

