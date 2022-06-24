function [ yRddot ] = yRefddot( t )
a=1;
b=0.25;
yRddot=-a*b^2*sin(b*t);


end

