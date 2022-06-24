function [ Dh_F ] = Dh_F_AUV_kinematic( U,X0,Udot,Xdot,h,P,N,Q,R,Qf,T,Umax,tau )

Dh_F=(1/h)*(F_AUV_kinematic( X0+h*Xdot,U+h*Udot,P,N,Q,R,Qf,T,Umax,tau )-F_AUV_kinematic( X0,U,P,N,Q,R,Qf,T,Umax,tau ) );

end

