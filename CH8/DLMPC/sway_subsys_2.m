function [ zeta_plus ] = sway_subsys_2( zeta,u,r,Fv,T )

x_wave=0;
y_wave=0;

% x_wave=0.2;
% y_wave=0.2;

rot_wave=0;

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

x = zeta(1);
y = zeta(2);
psi = zeta(3);
v = zeta(4);

u_wave = x_wave*cos(psi) + y_wave*sin(psi);
v_wave = -x_wave*sin(psi) + y_wave*cos(psi);

x_dot=(u+u_wave)*cos(psi)-(v+v_wave)*sin(psi);
y_dot=(u+u_wave)*sin(psi)+(v+v_wave)*cos(psi);

% x_dot = u*cos(psi)-v*sin(psi);
% y_dot = u*sin(psi)+v*cos(psi);
psi_dot = r;
v_dot=-(Mx/My)*u*r-(Yv/My)*v-(Dv/My)*v*abs(v)+Fv/My;

zeta_dot = [x_dot;y_dot;psi_dot;v_dot];
zeta_plus = zeta + zeta_dot*T;



end

