function [ tau ] = BS_controller( X0,t,Lambda,Kp,Kd )

% System coefficients
%==========================================================================
    err_model=0;

    m=116;
    Iz=13.1;
    X_udot=-167.6*(1+err_model);
    Y_vdot=-477.2*(1+err_model);
    N_rdot=-15.9*(1+err_model);
    Xu=26.9*(1+err_model);
    Yv=35.8*(1+err_model);
    Nr=3.5*(1+err_model);
    Du=241.3*(1+err_model);
    Dv=503.8*(1+err_model);
    Dr=76.9*(1+err_model);

    Mx=m-X_udot;
    My=m-Y_vdot;
    Mpsi=Iz-N_rdot;
%==========================================================================
    M=diag([Mx;My;Mpsi]);
    PI=[0,-1,0;1,0,0;0,0,0];
    g=0;

    xR=xRef(t);
    xRdot=xRefdot(t);
    xRddot=xRefddot(t); 
    xRdddot=xRefdddot(t); 
    
    yR=yRef(t);
    yRdot=yRefdot(t);
    yRddot=yRefddot(t);
    yRdddot=yRefdddot(t); 
    
    psiR=psiRef(xRdot,yRdot);
    psiRdot=psiRefdot(xRdot,xRddot,yRdot,yRddot);
    psiRddot=psiRefddot(xRdot,xRddot,xRdddot,yRdot,yRddot,yRdddot);

    eta=X0(1:3,1);
    vel=X0(4:6,1);
    psi=X0(3,1);
    
    eta_d=[xR;yR;psiR];
    eta_d_dot=[xRdot;yRdot;psiRdot];
    eta_d_ddot=[xRddot;yRddot;psiRddot];
    
    tilde_eta=eta-eta_d;
    eta_r_dot=eta_d_dot-Lambda*tilde_eta;
    
    Rpsi=Rot(psi);
    v_r=Rpsi'*eta_r_dot;
    
    eta_dot=Rpsi*vel;
    S=eta_dot-eta_r_dot;
    
    psi_dot=eta_dot(3);
    v_r_dot=-psi_dot*Rpsi'*PI*eta_r_dot+Rpsi'*(eta_d_ddot-Lambda*(eta_dot-eta_d_dot));
    
    C=CV(vel);
    D=DV(vel);
    
    tau=M*v_r_dot+C*v_r+D*v_r+g-Rpsi'*Kp*tilde_eta-Rpsi'*Kd*S;

end

