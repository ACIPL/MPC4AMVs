function [ Xplus ] = AUV_Dynamics_discrete( X,U,T )

x_wave=0;
y_wave=0;

% x_wave=0.2;
% y_wave=0.2;

rot_wave=0;

% System coefficients
%==========================================================================
err_model=0;

m=116;
Iz=13.1;
X_udot=-167.6*(1+err_model);
Y_vdot=-477.2*(1+err_model);
N_rdot=-15.9*(1+err_model);
Xu=26.9*(1+err_model);
Yv=35.8*(1+err_model);
Nr=3.5*(1+err_model);
Du=241.3*(1+err_model);
Dv=503.8*(1+err_model);
Dr=76.9*(1+err_model);

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

Fu=U(1);
Fv=U(2);
Fr=U(3);

u_wave = x_wave*cos(psi) + y_wave*sin(psi);
v_wave = -x_wave*sin(psi) + y_wave*cos(psi);

x_dot=(u+u_wave)*cos(psi)-(v+v_wave)*sin(psi);
y_dot=(u+u_wave)*sin(psi)+(v+v_wave)*cos(psi);

psi_dot=r;
u_dot=(My/Mx)*v*r-(Xu/Mx)*u-(Du/Mx)*u*abs(u)+Fu/Mx;
v_dot=-(Mx/My)*u*r-(Yv/My)*v-(Dv/My)*v*abs(v)+Fv/My;
r_dot=((Mx-My)/Mpsi)*u*v-(Nr/Mpsi)*r-(Dr/Mpsi)*r*abs(r)+Fr/Mpsi;


X_dot=[x_dot;y_dot;psi_dot;u_dot;v_dot;r_dot];

Xplus=X+X_dot*T;

end

