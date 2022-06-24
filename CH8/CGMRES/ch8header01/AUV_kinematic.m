classdef AUV_kinematic<handle
    
    properties
        model Model = AUV(ones(11,1),3,zeros(6,1),zeros(3,1))
        weights
    end
    
    methods
        function obj = AUV_kinematic(model,weights)
            obj.model = model;
            obj.weights = weights; % Q R Qf
        end
        
        function Hu = Hu_AUV_kinematic( obj, lambda,U,tau,Umax )
            % U=[Fu;Fv;Fr;u4;u5;u6;mu1;mu2;mu3]
            
            % System coefficients
            %==============================================================
            Mx = obj.model.DeducedCoef(1);
            My = obj.model.DeducedCoef(2);
            Mpsi = obj.model.DeducedCoef(3);
            %==============================================================
            R = obj.weights{2};
            r11=R(1,1);
            r22=R(2,2);
            r33=R(3,3);
            
            lambda1=lambda(1);
            lambda2=lambda(2);
            lambda3=lambda(3);
            lambda4=lambda(4);
            lambda5=lambda(5);
            lambda6=lambda(6);
            
            Fu=U(1);
            Fv=U(2);
            Fr=U(3);
            
            Fu_max=Umax(1);
            Fv_max=Umax(2);
            Fr_max=Umax(3);
            
            Hu1=r11*Fu+lambda4/Mx+tau/(Fu_max-Fu)-tau/(Fu_max+Fu);
            Hu2=r22*Fv+lambda5/My+tau/(Fv_max-Fv)-tau/(Fv_max+Fv);
            Hu3=r33*Fr+lambda6/Mpsi+tau/(Fr_max-Fr)-tau/(Fr_max+Fr);
            
            Hu=[Hu1;Hu2;Hu3];
        end
        %
        function Hx = Hx_AUV_kinematic(obj, X,lambda,U,p )
            
            % System coefficients
            %==============================================================
            Mx = obj.model.DeducedCoef(1);
            My = obj.model.DeducedCoef(2);
            Mpsi = obj.model.DeducedCoef(3);
            Xu = obj.model.Coef(6);
            Yv = obj.model.Coef(7);
            Nr = obj.model.Coef(8);
            Du = obj.model.Coef(9);
            Dv = obj.model.Coef(10);
            Dr = obj.model.Coef(11);
            %==============================================================
            
            x=X(1); y=X(2); psi=X(3);
            u=X(4); v=X(5); r=X(6);
            
            xd=p(1); yd=p(2); psi_d=p(3);
            ud=p(4); vd=p(5); rd=p(6);
            
            lambda1=lambda(1); lambda2=lambda(2); lambda3=lambda(3);
            lambda4=lambda(4); lambda5=lambda(5); lambda6=lambda(6);
            
            Fu=U(1);
            Fv=U(2);
            Fr=U(3);
            
            Q = obj.weights{1};
            q11=Q(1,1); q22=Q(2,2); q33=Q(3,3);
            q44=Q(4,4); q55=Q(5,5); q66=Q(6,6);
            
            Hx1=q11*(x-xd);
            Hx2=q22*(y-yd);
            Hx3=q33*(psi-psi_d)-lambda1*u*sin(psi)-lambda1*v*cos(psi)+lambda2*u*cos(psi)-lambda2*v*sin(psi);
            Hx4=q44*(u-ud)+lambda1*cos(psi)+lambda2*sin(psi)-lambda4*(Xu/Mx)-lambda5*(Mx/My)*r+lambda6*((Mx-My)/Mpsi)*v-2*lambda4*(Du/Mx)*abs(u);
            Hx5=q55*(v-vd)-lambda1*sin(psi)+lambda2*cos(psi)+lambda4*(My/Mx)*r-lambda5*(Yv/My)+lambda6*((Mx-My)/Mpsi)*u-2*lambda5*(Dv/My)*abs(v);
            Hx6=q66*(r-rd)+lambda3+lambda4*(My/Mx)*v-lambda5*(Mx/My)*u-lambda6*(Nr/Mpsi)-2*lambda6*(Dr/Mpsi)*abs(r);

            Hx=[Hx1;Hx2;Hx3;Hx4;Hx5;Hx6];
            
        end
        
        function PHIx = PHIx_AUV_kinematic(obj, X, pf )
            x=X(1); y=X(2); psi=X(3);
            u=X(4); v=X(5); r=X(6);
            
            xf=pf(1);yf=pf(2); psi_f=pf(3);
            uf=pf(4); vf=pf(5); rf=pf(6);
            
            Qf = obj.weights{3};
            qf11=Qf(1,1); qf22=Qf(2,2); qf33=Qf(3,3);
            qf44=Qf(4,4); qf55=Qf(5,5); qf66=Qf(6,6);
            
            
            PHIx1=qf11*(x-xf);
            PHIx2=qf22*(y-yf);
            PHIx3=qf33*(psi-psi_f);
            PHIx4=qf44*(u-uf);
            PHIx5=qf55*(v-vf);
            PHIx6=qf66*(r-rf);
            
            PHIx=[PHIx1;PHIx2;PHIx3;PHIx4;PHIx5;PHIx6];
            
        end
        
        function F = F_AUV_kinematic(obj, X0,U,P,N,T,Umax,tau )
            % U has size of [9,N]; U=[U0,U1,...,U_N-1]
            % Ui=[Fu;Fv;Fr;u4;u5;u6;mu1;mu2;mu3]
            
            % Umax=[Fu_max;Fv_max;Fr_max]
            
            % P has size of [6,N+1]; P=[Xd0,Xd1,...,Xd_N]
            % Xdi=[xd;yd;psi_d;ud;vd;rd]
            
            % X evolution throuth system dynamics
            %==============================================================
            Q = obj.weights{1};
            R = obj.weights{2};
            Qf = obj.weights{3};
            
            [nU,mU]=size(U);
            
            m=length(X0);
            X=zeros(m,N+1); % X=[X0,X1,...,XN]
            
            X(:,1)=X0;
            for j=1:1:N
                X(:,j+1)=obj.model.dynamics_discrete( X(:,j),U(:,j),T);
            end
            %==============================================================
            %Lambda propagation backwardly
            %==============================================================
            Lambda=zeros(m,N+1); % Lambda=[lambda0,lambda1,...,lambdaN]
            Lambda(:,N+1)=obj.PHIx_AUV_kinematic( X(:,N+1),P(:,N+1));
            
            for i=N:-1:1
                Hx = obj.Hx_AUV_kinematic( X(:,i),Lambda(:,i+1),U,P(:,i));
                Lambda(:,i)=Lambda(:,i+1)+Hx*T;
            end
            %==============================================================
            F=zeros(nU*N,1);
            
            for k=1:1:N
                F1=obj.Hu_AUV_kinematic( Lambda(:,k+1),U(:,k),tau,Umax);
                F(nU*(k-1)+1:nU*k,1)=F1;
            end
        end
        function U_dot = FDGMRES_AUV_kinematic(obj, U,X0,Xdot,h,P,N,T,U_dot_hat,k_max,tol,As,Umax,tau )
            
            [nU,mU]=size(U);
            m=length(X0);
            
            b=As*obj.F_AUV_kinematic( X0,U,P,N,T,Umax,tau )-obj.Dh_F_AUV_kinematic( U,X0,zeros(nU,N),Xdot,h,P,N,T,Umax,tau );
            r_hat=b-obj.Dh_F_AUV_kinematic( U,X0+h*Xdot,U_dot_hat,zeros(m,1),h,P,N,T,Umax,tau );
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
                Vk(:,k+1)= obj.Dh_F_AUV_kinematic( U,X0+h*Xdot,Vkk_reshape,zeros(m,1),h,P,N,T,Umax,tau );
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
            k_out;
            yk;
            Vk;
            U_dot_hat_reshape=reshape(U_dot_hat,[nU*N,1]);
            U_dot=U_dot_hat_reshape+Vk(:,1:k_out)*yk;
            U_dot=reshape(U_dot,[nU,N]);
        end
        
        function Dh_F = Dh_F_AUV_kinematic( obj,U,X0,Udot,Xdot,h,P,N,T,Umax,tau )
            Dh_F=(1/h)*(obj.F_AUV_kinematic( X0+h*Xdot,U+h*Udot,P,N,T,...
                Umax,tau )-obj.F_AUV_kinematic( X0,U,P,N,T,Umax,tau ) );
        end
        
        
    end
end

