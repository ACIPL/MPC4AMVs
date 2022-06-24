function [ Hu ] = Hu_AUV_kinematic( X,lambda,U,R,tau,Umax )
% U=[Fu;Fv;Fr;u4;u5;u6;mu1;mu2;mu3]

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

r11=R(1,1);
r22=R(2,2);
r33=R(3,3);

lambda1=lambda(1);
lambda2=lambda(2);
lambda3=lambda(3);
lambda4=lambda(4);
lambda5=lambda(5);
lambda6=lambda(6);

Fu=U(1);
Fv=U(2);
Fr=U(3);

Fu_max=Umax(1);
Fv_max=Umax(2);
Fr_max=Umax(3);

Hu1=r11*Fu+lambda4/Mx+tau/(Fu_max-Fu)-tau/(Fu_max+Fu);
Hu2=r22*Fv+lambda5/My+tau/(Fv_max-Fv)-tau/(Fv_max+Fv);
Hu3=r33*Fr+lambda6/Mpsi+tau/(Fr_max-Fr)-tau/(Fr_max+Fr);

Hu=[Hu1;Hu2;Hu3];

end

