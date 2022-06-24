function [ C ] = CV( vel )

% System coefficients
%==========================================================================
m = 116;
X_udot = -167.6;
Y_vdot = -477.2;
Mudot = m - X_udot;
Mvdot = m - Y_vdot;
%==========================================================================
u = vel(1);
v = vel(2);

C1=[0;0;Mvdot*v];
C2=[0;0;-Mudot*u];
C3=[-Mvdot*v;Mudot*u;0];

C=[C1,C2,C3];

end

