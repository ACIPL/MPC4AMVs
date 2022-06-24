classdef MPC_controller < handle
    properties
        N
        model Model = AUVI(ones(11,1),3,zeros(6,1),zeros(3,1))
        weights
        upperbound
        lowerbound
    end
    
    methods
        function obj = MPC_controller(N,model,weights,upperbound,lowerbound)
            obj.N = N;
            obj.model = model;
            obj.weights = weights;
            obj.calc_upperbound(upperbound);
            obj.calc_lowerbound(lowerbound);
        end
        
        function obj = calc_upperbound(obj,upperbound)
            if ~isempty(upperbound)
                m = length(upperbound);
                Thrust_max = zeros(m*obj.N,1);
                for i=1:1:obj.N
                    Thrust_max(1+m*(i-1):m*i,1) = upperbound;
                end
                obj.upperbound = Thrust_max;
            else
                obj.upperbound = [];
            end
        end
        
        function obj = calc_lowerbound(obj,lowerbound)
            if ~isempty(lowerbound)
                m = length(lowerbound);
                Thrust_min = zeros(m*obj.N,1);
                for i=1:1:obj.N
                    Thrust_min(1+m*(i-1):m*i,1) = lowerbound;
                end
                obj.lowerbound = Thrust_min;
            else
                obj.lowerbound = [];
            end
        end
        
        function cost = mpc_cost(obj,u,X0,zeta0,dt)
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
            Hp = obj.N;
            Q = obj.weights{1}(1:6,1:6);
            R = obj.weights{2}(1:3,1:3);
            Qf = obj.weights{3}(1:6,1:6);
            
            nu = length(u);
            nu = nu/Hp;
            nx = length(X0);
            Hu = Hp;
            U = zeros(nu,Hp);
            X = zeros(nx,Hp);
            ZETA=zeros(2,Hp); % zeta=[s,s_dot];
            
            xr = zeros(Hp,1);
            yr = zeros(Hp,1);
            psi_r = zeros(Hp,1);
            ur = zeros(Hp,1);
            vr = zeros(Hp,1);
            rr = zeros(Hp,1);
            P=zeros(9,Hp); % P(:,i) = [xri;yri;psi_ri;uri;vri;rri;urdi;vrdi;rrdi];
            
            ur_d=zeros(Hp,1);
            rr_d=zeros(Hp,1);
            
            Fur = zeros(Hp,1);
            Fvr = zeros(Hp,1);
            Frr = zeros(Hp,1);
            vsr=zeros(Hp,1);
            
            %==========================================================================%
            % partition of u
            %==========================================================================%
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i) = u((i-1)*nu+j,1);
                end
            end
            
            %==========================================================================%
            % state prediction
            %==========================================================================%
            Xplus = obj.model.dynamics_discrete( X0,U(:,1),dt );
            X(:,1) = Xplus;
            zeta_plus = path_evol_2o(zeta0,U(4,1),dt);
            ZETA(:,1)=zeta_plus;
            
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,U(1:3,i),dt );
                X(:,i)= Xplus;
                zeta_plus = path_evol_2o( zeta_plus,U(4,i),dt );
                ZETA(:,i) = zeta_plus;
            end
            
            %==========================================================================%
            % reference states
            %==========================================================================%
            % reference state 1
            for i=1:1:Hp
                xr(i,1)=ZETA(1,i);
                yr(i,1)=sin(ZETA(1,i));
                psi_r(i,1)=atan(cos(ZETA(1,i)));
                
                ur(i,1)=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)^(1/2)*ZETA(2,i);
                vr(i,1)=0;
                rr(i,1)=ZETA(2,i)*(px_dot( ZETA(1,i) )*py_ddot( ZETA(1,i) )...
                    -py_dot( ZETA(1,i) )*px_ddot( ZETA(1,i) ))/...
                    (px_dot( ZETA(1,i) )^2+py_dot( ZETA(1,i) )^2);
                
                vr_d = 0;
                
                urd1=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)^(-1/2)*(px_dot(ZETA(1,i))*px_ddot(ZETA(1,i))+py_dot(ZETA(1,i))*py_ddot(ZETA(1,i)))*U(4,i);
                urd2=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)^(1/2)*U(4,i);
                ur_d(i,1)=urd1+urd2;
                
                rrd1=(px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2)*(px_dot(ZETA(1,i))*py_dddot(ZETA(1,i))-py_dot(ZETA(1,i))*px_dddot(ZETA(1,i)))-2*(px_dot(ZETA(1,i))*py_ddot(ZETA(1,i))-py_dot(ZETA(1,i))*px_ddot(ZETA(1,i)))*(px_dot(ZETA(1,i))*px_ddot(ZETA(1,i))+py_dot(ZETA(1,i))*py_ddot(ZETA(1,i)));
                rrd2=px_dot(ZETA(1,i))*py_ddot(ZETA(1,i))-py_dot(ZETA(1,i))*px_ddot(ZETA(1,i));
                deno=px_dot(ZETA(1,i))^2+py_dot(ZETA(1,i))^2;
                rr_d(i,1)=(rrd1/deno^2+rrd2/deno)*U(4,i);
                
                P(:,i) = [xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1);ur_d(i,1);vr_d;rr_d(i,1)];
                
                Fur(i,1)=Mx*ur_d(i,1)+Xu*ur(i,1)+Du*ur(i,1)*abs(ur(i,1));
                Fvr(i,1)= Mx*ur(i,1)*rr(i,1);
                Frr(i,1)=Mpsi*rr_d(i,1)+Nr*rr(i,1)+Dr*rr(i,1)*abs(rr(i,1));
                vsr(i,1)=U(4,i);
            end
            
            
            %==========================================================================%
            % cost function
            %==========================================================================%
            cost = 0;
            
            for i=1:1:Hp
                
                Xr=[xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1)];
                Xe = ErrorState( X(:,i), Xr );
                f1 = F1( X(:,i),Xe,P(:,i) );
                f2 = F2( X(:,i),Xe,P(:,i) );
                f3 = F3( X(:,i),Xe,P(:,i) );
                tau_u = (U(1,i)-f1)/Mx;
                tau_v = (U(2,i)-f2)/My;
                tau_r = (U(3,i)-f3)/Mpsi;
                Te = [tau_u;tau_v;tau_r];
                cost = cost+Xe'*Q*Xe+Te'*R*Te;
                
                %      fx = fx+q77*(ZETA(1,i)-ZETA(1,i))^2 + r44*(U(4,i)-vsr(i,1))^2; % equal to zero
            end
            
            Xr=[xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
            Xe = ErrorState( X(:,Hp), Xr );
            cost = cost+Xe'*Qf*Xe;
        end
        function cost = mpc_cost2(obj,u,X0,zeta0,dt)
            
            s0 = zeta0(1,1);
            Hp = obj.N;
            Hu = Hp;
            nu=length(u);
            nu=nu/Hp;
            nx=length(X0);
            
            Q = obj.weights{4};
            R = obj.weights{5};
            Qf = obj.weights{6};
            U=zeros(nu,Hp);
            X=zeros(nx,Hp);
            ZETA=zeros(2,Hp);
            S=zeros(1,Hp);
            %==========================================================================
            % partition of u
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i)=u((i-1)*nu+j,1);
                end
            end
            %==========================================================================
            % state prediction
            Xplus = obj.model.dynamics_discrete( X0,U(1:3,1),dt );
            X(:,1)=Xplus;
            zeta_plus = path_evol_2o( zeta0,U(4,1),dt );
            ZETA(:,1)=zeta_plus;
            S(1)=ZETA(1,1);
            
            for i=2:1:Hp
                Xplus =obj.model.dynamics_discrete( Xplus,U(1:3,i),dt );
                X(:,i)= Xplus;
                
                zeta_plus = path_evol_2o( zeta_plus,U(4,i),dt );
                ZETA(:,i) = zeta_plus;
                S(i) = ZETA(1,i);
            end
            %==========================================================================
            % reference state generation
            [ Xr_N,Sr_N,Ur_N ] = Ref_State_Mat( s0,Hp,dt );
            %==========================================================================
            % cost function
            %uR=1;
            cost=0;
            for i=1:1:Hp
                Xaug=[X(:,i);S(i)];
                Xaug_r=[Xr_N(:,i);Sr_N(i)];
                cost=cost+(Xaug_r-Xaug)'*Q*(Xaug_r-Xaug);     
            end
            
            for i=1:1:Hu
                cost=cost+(Ur_N(:,i)-U(:,i))'*R*(Ur_N(:,i)-U(:,i));     
            end
            
            Xaugf=[X(:,Hp);S(Hp)];
            Xaugf_r=[Xr_N(:,Hp);Sr_N(Hp)];
            cost=cost+(Xaugf_r-Xaugf)'*Qf*(Xaugf_r-Xaugf);
        end
        
        function [y,ceq] = mpc_contraints(obj,u,X0,zeta0,dt)
            nx=length(X0);
            nu=length(u);
            Hp = obj.N;
            nu=nu/Hp;
            U=zeros(nu,Hp);
            
            %==========================================================================
            % partition of u 
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i)=u((i-1)*nu+j,1);
                end
            end
            
            %==========================================================================
            %state prediction
            X=zeros(nx,Hp);
            Xplus=X0;
            
            ZETA=zeros(2,Hp); % zeta=[s,s_dot];
            zeta_plus=zeta0;
            y=zeros(Hp,1);
            
            for i=1:1:Hp
                Xplus =obj.model.dynamics_discrete( Xplus,U(1:3,i),dt );
                X(:,i)= Xplus;
                zeta_plus = path_evol_2o( zeta_plus,U(4,i),dt );
                ZETA(:,i)=zeta_plus;
                y(i)=0-ZETA(2,i);
            end
            
            xrN=ZETA(1,Hp);
            yrN=sin(ZETA(1,Hp));
            psi_rN=atan(cos(ZETA(1,Hp)));
            
            urN=(px_dot(ZETA(1,Hp))^2+py_dot(ZETA(1,Hp))^2)^(1/2)*ZETA(2,Hp);
            vrN=0;
            rrN=ZETA(2,Hp)*(px_dot( ZETA(1,Hp) )*py_ddot( ZETA(1,Hp) )-py_dot( ZETA(1,Hp) )*px_ddot( ZETA(1,Hp) ))/(px_dot( ZETA(1,Hp) )^2+py_dot( ZETA(1,Hp) )^2);
            
            XrN=[xrN,yrN,psi_rN,urN,vrN,rrN]';
            ceq=X(:,Hp)-XrN;
        end
        
        function [y,ceq] = mpc_contraints2(obj,u,J1_star,epsi,X0,zeta0,dt)
            nx=length(X0);
            nu=length(u);
            Hp = obj.N;
            nu=nu/Hp;
            U=zeros(nu,Hp);
            
            %==========================================================================
            % partition of u
            
            for i=1:1:Hp
                for j=1:1:nu
                    
                    U(j,i)=u((i-1)*nu+j,1);
                    
                end
            end
            
            %==========================================================================
            %state prediction
            
            X=zeros(nx,Hp);
            Xplus=X0;
            
            ZETA=zeros(2,Hp); % zeta=[s,s_dot];
            zeta_plus=zeta0;
            y1=zeros(Hp,1);
            
            for i=1:1:Hp
                
                Xplus =obj.model.dynamics_discrete( Xplus,U(1:3,i),dt );
                X(:,i)= Xplus;
                zeta_plus = path_evol_2o( zeta_plus,U(4,i),dt );
                ZETA(:,i)=zeta_plus;
                y1(i)=0-ZETA(2,i);
                
                
            end
            
            xrN = ZETA(1,Hp);
            yrN = sin(ZETA(1,Hp));
            psi_rN = atan(cos(ZETA(1,Hp)));
            
            urN = (px_dot(ZETA(1,Hp))^2+py_dot(ZETA(1,Hp))^2)^(1/2)*ZETA(2,Hp);
            vrN = 0;
            rrN = ZETA(2,Hp)*(px_dot( ZETA(1,Hp) )*py_ddot( ZETA(1,Hp) )-py_dot( ZETA(1,Hp) )*px_ddot( ZETA(1,Hp) ))/(px_dot( ZETA(1,Hp) )^2+py_dot( ZETA(1,Hp) )^2);
            
            XrN = [xrN,yrN,psi_rN,urN,vrN,rrN]';
            
            ceq = X(:,Hp)-XrN;
            
            J1 = obj.mpc_cost(u,X0,zeta0,dt);
            y2 = J1-J1_star-epsi;
            
            y = [y1;y2];
            
        end
        function [u, J1] = calc_control1(obj,zeta0,X0,u0,dt)
            options = optimset('Algorithm','sqp');
            [u, J1] = fmincon(@(u) obj.mpc_cost(u,X0,zeta0,dt),u0,[],[],...
                [],[],obj.lowerbound,obj.upperbound,...
                @(u)obj.mpc_contraints(u,X0,zeta0,dt),options);
        end 
        function [u, J2] = calc_control2(obj,J1_star,epsi,zeta0,X0,u0,dt)
            options = optimset('Algorithm','sqp');
            [u, J2] = fmincon(@(u) obj.mpc_cost2(u,X0,zeta0,dt),u0,[],...
                [],[],[],obj.lowerbound,obj.upperbound,...
                @(u)obj.mpc_contraints2(u,J1_star,epsi,X0,zeta0,dt),options);
        end 
        
        %%%%%%%%%%%%%%%%%%
        
        function get_params(obj)
            fprintf('weights: \n')
            disp(obj.weights)
        end
    end
end


