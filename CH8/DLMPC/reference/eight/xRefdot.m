function [ xRdot ] = xRefdot( t )
a=1;
b=-0.5;

xRdot=a*b*cos(b*t);


end

