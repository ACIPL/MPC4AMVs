function [ Xe ] = ErrorState( X, Xr )

x = X(1);
y = X(2);
psi = X(3);
u = X(4);
v = X(5);
r = X(6);

xr = Xr(1);
yr = Xr(2);
psi_r = Xr(3);
ur = Xr(4);
vr = Xr(5);
rr = Xr(6);

xe = (xr-x)*cos(psi)+(yr-y)*sin(psi);
ye = -(xr-x)*sin(psi)+(yr-y)*cos(psi);
psi_e = psi_r-psi;

ue = u-ur*cos(psi_e)+vr*sin(psi_e);
ve = v-ur*sin(psi_e)-vr*cos(psi_e);
re = r-rr;

Xe = [xe;ye;psi_e;ue;ve;re];
end

