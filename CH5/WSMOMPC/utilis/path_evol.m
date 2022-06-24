function [ s_plus ] = path_evol( s,v_s,T )
lambda=0;
s_dot=-lambda*s+v_s;
s_plus=s+s_dot*T;

end

