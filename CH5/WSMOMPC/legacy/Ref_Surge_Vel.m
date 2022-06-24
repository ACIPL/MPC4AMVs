function [ v_sr ] = Ref_Surge_Vel( s, uR )

 v_sr=px_dot( s )^2+py_dot( s )^2;
 
 v_sr=v_sr^(-1/2)*uR;


end

