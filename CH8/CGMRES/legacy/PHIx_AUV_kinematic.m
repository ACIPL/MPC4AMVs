function [ PHIx ] = PHIx_AUV_kinematic( X,Qf,pf )
x=X(1);
y=X(2);
psi=X(3);
u=X(4);
v=X(5);
r=X(6);

xf=pf(1);
yf=pf(2);
psi_f=pf(3);
uf=pf(4);
vf=pf(5);
rf=pf(6);


qf11=Qf(1,1);
qf22=Qf(2,2);
qf33=Qf(3,3);
qf44=Qf(4,4);
qf55=Qf(5,5);
qf66=Qf(6,6);


PHIx1=qf11*(x-xf);
PHIx2=qf22*(y-yf);
PHIx3=qf33*(psi-psi_f);
PHIx4=qf44*(u-uf);
PHIx5=qf55*(v-vf);
PHIx6=qf66*(r-rf);

PHIx=[PHIx1;PHIx2;PHIx3;PHIx4;PHIx5;PHIx6];

end

