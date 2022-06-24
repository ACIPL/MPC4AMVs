function [ zeta_plus ] = path_evol_2o( zeta,v_s,T )
% zeta=[s s_dot]';

A=[0,1;0,0];
B=[0;1];

zeta_d=A*zeta+B*v_s;
zeta_plus=zeta+zeta_d*T;

end

