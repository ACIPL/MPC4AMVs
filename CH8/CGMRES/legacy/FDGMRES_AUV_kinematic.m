function [ U_dot ] = FDGMRES_AUV_kinematic( U,X0,Xdot,h,P,N,Q,R,Qf,T,U_dot_hat,k_max,tol,As,Umax,tau )

[nU,mU]=size(U);
m=length(X0);

b=As*F_AUV_kinematic( X0,U,P,N,Q,R,Qf,T,Umax,tau )-Dh_F_AUV_kinematic( U,X0,zeros(nU,N),Xdot,h,P,N,Q,R,Qf,T,Umax,tau );
r_hat=b-Dh_F_AUV_kinematic( U,X0+h*Xdot,U_dot_hat,zeros(m,1),h,P,N,Q,R,Qf,T,Umax,tau );
v1=r_hat/norm(r_hat);
nv=length(v1);
Vk=zeros(nv,k_max);
rho=norm(r_hat);
beta=rho;
k=0;
Vk(:,1)=v1;

Hk=zeros(k_max+1,k_max);

while k<k_max
    k=k+1;
    Vkk_reshape=reshape(Vk(:,k),[nU,N]);
    Vk(:,k+1)=Dh_F_AUV_kinematic( U,X0+h*Xdot,Vkk_reshape,zeros(m,1),h,P,N,Q,R,Qf,T,Umax,tau );
    for j=1:1:k
        Hk(j,k)=Vk(:,k+1)'*Vk(:,j);
        Vk(:,k+1)=Vk(:,k+1)-Hk(j,k)*Vk(:,j);
    end
    Hk(k+1,k)=norm(Vk(:,k+1));
    Vk(:,k+1)=Vk(:,k+1)/norm(Vk(:,k+1));
    e1=[1;zeros(k,1)];
    Hk_hat=Hk(1:k+1,1:k);
    
    b_hat=beta*e1;
    yk=pinv(Hk_hat)*b_hat;
    rho=norm(b_hat-Hk_hat*yk);
    k_out=k;
    if rho<=tol
       break;
    end
    
end
k_out
yk;
Vk;
U_dot_hat_reshape=reshape(U_dot_hat,[nU*N,1]);
U_dot=U_dot_hat_reshape+Vk(:,1:k_out)*yk;

U_dot=reshape(U_dot,[nU,N]);

end

