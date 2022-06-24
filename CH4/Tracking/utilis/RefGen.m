function [ P ] = RefGen(M,N,t0,dt,coef,Bplus,traj)
    m = coef(1);
    Iz = coef(2);
    X_udot = coef(3);
    N_rdot = coef(5);
    Xu = coef(6);
    Nr = coef(8);
    Du = coef(9);
    Dr = coef(11);
    Mudot = m - X_udot;
    Mrdot = Iz - N_rdot;
  
    P = zeros(M,N);
    
    for j=1:1:N % Generate reference upto N steps ahead
        xR = traj.Xr(t0);
        xRdot = traj.Xrdot(t0);
        xRddot = traj.Xrddot(t0);
        yR = traj.Yr(t0);
        yRdot = traj.Yrdot(t0);
        yRddot = traj.Yrddot(t0);
        psiR = traj.PSIr(t0);
        psiRdot = traj.PSIrdot(t0);
        psiRddot = traj.PSIrddot(t0);      
        uR = sqrt(xRdot^2+yRdot^2);
        vR = 0;
        rR = psiRdot;   
        uRdot = (xRdot^2+yRdot^2)^(-1/2)*(xRdot*xRddot+yRdot*yRddot);
        rRdot = psiRddot;        
        FuR = Mudot*uRdot+Xu*uR+Du*uR*abs(uR);
        FvR = Mudot*uR*rR;
        FrR = Mrdot*rRdot+Nr*rR+Dr*rR*abs(rR);
        Tau = [FuR;FvR;FrR];
        Thrust = Bplus*Tau;
        P(:,j)=[xR;yR;psiR;uR;vR;rR;FuR;FvR;FrR;Thrust];
        t0 = t0 + dt;
    end

end

