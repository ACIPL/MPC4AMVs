function [ C ] = CV( vel )

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

u=vel(1);
v=vel(2);
r=vel(3);

C1=[0;0;My*v];
C2=[0;0;-Mx*u];
C3=[-My*v;Mx*u;0];

C=[C1,C2,C3];

end

