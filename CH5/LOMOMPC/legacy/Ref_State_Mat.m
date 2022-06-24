function [ Xr_N,Sr_N,Ur_N ] = Ref_State_Mat( s0,uR,N,T )

% AUV System coefficients
%==========================================================================
m=116;
Iz=13.1;
X_udot=-167.6;
Y_vdot=-477.2;
N_rdot=-15.9;
Xu=26.9;
Yv=35.8;
Nr=3.5;
Du=241.3;
Dv=503.8;
Dr=76.9;
Mx=m-X_udot;
My=m-Y_vdot;
Mpsi=Iz-N_rdot;
%==========================================================================

Xr_N=zeros(6,N);
Sr_N=zeros(N,1);
Ur_N=zeros(4,N);

v_sr = Ref_Surge_Vel( s0, uR );
sr_plus=path_evol( s0,v_sr,T );

for i=1:1:N
    xr=px( sr_plus );
    yr=py( sr_plus );
    psi_r=atan2(py_dot( sr_plus ),px_dot( sr_plus ));
    u_r=uR;
    v_r=0;
    r_r=v_sr*(px_dot( sr_plus )*py_ddot( sr_plus )-py_dot( sr_plus )*px_ddot( sr_plus ))/(px_dot( sr_plus )^2+py_dot( sr_plus )^2);
    Xr=[xr;yr;psi_r;u_r;v_r;r_r];
    Xr_N(:,i)=Xr;
    
    Sr_N(i)=sr_plus;
       
    Fur=-My*v_r*r_r+Xu*u_r+Du*u_r*abs(u_r);
    Fvr= Mx*u_r*r_r+Yv*v_r+Dv*v_r*abs(v_r);
    Frr=(My-Mx)*u_r*v_r+Nr*r_r+Dr*r_r*abs(r_r)+ u_r*(-cos(sr_plus)*(1+cos(sr_plus)^2)^(3/2)+3*(1+cos(sr_plus)^2)^(1/2)*sin(sr_plus)^2*cos(sr_plus))/3*(1+cos(sr_plus)^2)^3;
    Ur=[Fur;Fvr;Frr;v_sr];
    Ur_N(:,i)=Ur;
    
    v_sr = Ref_Surge_Vel( sr_plus, uR );
    sr_plus=path_evol( sr_plus,v_sr,T );
    
end



end

