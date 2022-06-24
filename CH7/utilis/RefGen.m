function P = RefGen(nx,N,t0,dt,X0)
P = zeros(nx,N);
Xplus = Ref( X0,dt, t0 );
P(:,1) = Xplus;
t0 = t0+dt;
for i=2:1:N
    Xplus = Ref( Xplus,dt, t0);
    P(:,i)= Xplus;
    t0 = t0 + dt;
end
end

function xd = Ref(xd, dt, t)
if t <= 2 && t > 1.6
    ad = 1.4;
elseif t > 3.9 && t <= 4.3
    ad = -1.6;
else
    ad = 0;
end

dxd = zeros(2,1);
dxd(1) = xd(2);
dxd(2) = ad;
xd = xd + dxd * dt;
end