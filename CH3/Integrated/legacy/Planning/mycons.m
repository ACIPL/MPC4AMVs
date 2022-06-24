function [ y,ceq ] = mycons( a,bbrks,mlti,T,X,Y,stroke )

a=a';
knots=augknt(bbrks,mlti);
sp=spmak(knots,a);

nbrk=length(bbrks);

x=bbrks(1):T:bbrks(nbrk);
x=x';
nx=length(x);
y=zeros(2*nx,1);
ceq=zeros(2*nx,1);

for i=1:1:nx;

    y(i)=fnval(sp,x(i))-plygfun2(x(i),X,Y,stroke);

end

for j=(nx+1):1:(2*nx)

    y(j)=plygfun1(x(j-nx),X,Y)-fnval(sp,x(j-nx));
end

end

