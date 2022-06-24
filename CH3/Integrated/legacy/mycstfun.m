function [ fx ] = mycstfun( a,bbrks,mlti,T )
a=a';
knots=augknt(bbrks,mlti);
sp=spmak(knots,a);

fprime=fnder(sp);
f2prime=fnder(sp,2);
f3prime=fnder(sp,3);

nbrk=length(bbrks);

fx=0;

for x=bbrks(1):T:bbrks(nbrk)

kprimex=(fnval(f3prime,x)*(1+fnval(fprime,x)^2)^(3/2)-fnval(f2prime,x)*(3/2)*(1+fnval(fprime,x)^2)^(1/2)*2*fnval(fprime,x)*fnval(f2prime,x))/(1+fnval(fprime,x)^2)^3;
kprimex=kprimex^2;
f0x=kprimex/(1+fnval(fprime,x)^2)^(1/2);
fx=fx+T*f0x;

end

end

