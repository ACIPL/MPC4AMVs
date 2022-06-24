function [ y,ceq ] = mycons_continuity( a,bbrks,mlti,T,X,Y,stroke,y0,y0d,y02d,eps )

a=a';
knots=augknt(bbrks,mlti);
sp=spmak(knots,a);

nbrk=length(bbrks);

x=bbrks(1):T:bbrks(nbrk);
x=x';
nx=length(x);
y=zeros(2*nx+6,1);
ceq=zeros(2*nx,1);

for i=1:1:nx;

    y(i)=fnval(sp,x(i))-plygfun2(x(i),X,Y,stroke);

end

for j=(nx+1):1:(2*nx)

    y(j)=plygfun1(x(j-nx),X,Y)-fnval(sp,x(j-nx));
end


y(2*nx+1)=fnval(sp,X(1))-y0-eps;
y(2*nx+2)=y0-fnval(sp,X(1))-eps;

y(2*nx+3)=fnval(fnder(sp,1),X(1))-y0d-eps;
y(2*nx+4)=y0d-fnval(fnder(sp,1),X(1))-eps;

y(2*nx+5)=fnval(fnder(sp,2),X(1))-y02d-eps;
y(2*nx+6)=y02d-fnval(fnder(sp,2),X(1))-eps;

end

