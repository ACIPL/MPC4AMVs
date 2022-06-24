function [ yRdot ] = yRefdot( t )

a=1;
b=0.25;

yRdot=a*b*cos(b*t);

end

