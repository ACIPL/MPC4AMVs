function [ xRdddot ] = xRefdddot( t )
a=1;
b=-0.5;

xRdddot=-a*b^3*cos(b*t);

end

