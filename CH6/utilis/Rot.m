function [ Rpsi ] = Rot( psi )

R1 = [cos(psi);sin(psi);0];
R2 = [-sin(psi);cos(psi);0];
R3 = [0;0;1];

Rpsi = [R1,R2,R3];

end

