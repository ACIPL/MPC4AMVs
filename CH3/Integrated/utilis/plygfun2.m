function [ y ] = plygfun2( x,X,Y,stroke )

Y=Y+stroke;
Big=find(X>=x);
Small=find(X<x);
nsml=length(Small);

if nsml==0
    y=Y(1);
else
Xup=X(Big(1));
Xlw=X(Small(nsml));
Yup=Y(Big(1));
Ylw=Y(Small(nsml));
y=(Yup-Ylw)/(Xup-Xlw)*(x-Xlw)+Ylw;
end

end

