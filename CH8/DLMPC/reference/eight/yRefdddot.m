function [ yRdddot ] = yRefdddot( t )
a=1;
b=0.25;

yRdddot=-a*b^3*cos(b*t);

end

