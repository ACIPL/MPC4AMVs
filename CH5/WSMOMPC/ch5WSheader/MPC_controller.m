classdef MPC_controller < Controller
    properties
        N
        model Model = AUVI(ones(11,1),3,zeros(6,1),zeros(3,1))
        weights
        upperbound
        lowerbound
        path
    end
    
    methods
        function obj = MPC_controller(N,model,weights,upperbound,lowerbound,path)
            obj.N = N;
            obj.model = model;
            obj.weights = weights;
            obj.calc_upperbound(upperbound);
            obj.calc_lowerbound(lowerbound);
            obj.path = path;
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
        
        function cost = mpc_cost(obj,uR,alpha,u,X0,s0,dt)
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
            Q1 = obj.weights{1};
            R1 = obj.weights{2};
            Qf1 = obj.weights{3};
            Q2 = obj.weights{4};
            R2 = obj.weights{5};
            Qf2 = obj.weights{6};
            
            nu = length(u);
            nu = nu/Hp;
            nx = length(X0);
            Hu = Hp;
            U = zeros(nu,Hp);
            X = zeros(nx,Hp);
            S = zeros(1,Hp);
            
            xr = zeros(Hp,1);
            yr = zeros(Hp,1);
            psi_r = zeros(Hp,1);
            ur = zeros(Hp,1);
            vr = zeros(Hp,1);
            rr = zeros(Hp,1);
            
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
            s_plus = path_evol(s0,U(4,1),dt);
            S(1)=s_plus;
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,U(:,i),dt );
                X(:,i)= Xplus;
                s_plus = path_evol( s_plus,U(4,i),dt );
                S(i)=s_plus;
            end
            
            %==========================================================================%
            % reference states
            %==========================================================================%
            % reference state 1
            for i=1:1:Hp
                xr(i,1)=S(i);
                yr(i,1)=sin(S(i));
                psi_r(i,1)=atan(cos(S(i)));
                
                ur(i,1)=(obj.path.px_dot(S(i))^2+obj.path.py_dot(S(i))^2)^(1/2)*U(4,i);
                vr(i,1)=0;
                rr(i,1)=U(4,i)*(obj.path.px_dot( S(i) )*obj.path.py_ddot( S(i) )...
                    -obj.path.py_dot( S(i) )*obj.path.px_ddot( S(i) ))/...
                    (obj.path.px_dot( S(i) )^2+obj.path.py_dot( S(i) )^2);
                
                Fur(i,1)=-My*vr(i,1)*rr(i,1)+Xu*ur(i,1)+Du*ur(i,1)*abs(ur(i,1));
                Fvr(i,1)= Mx*ur(i,1)*rr(i,1)+Yv*vr(i,1)+Dv*vr(i,1)*abs(vr(i,1));
                Frr(i,1)=(My-Mx)*ur(i,1)*vr(i,1)+Nr*rr(i,1)+Dr*rr(i,1)*abs(rr(i,1))+ ur(i,1)*(-cos(S(i))*(1+cos(S(i))^2)^(3/2)+3*(1+cos(S(i))^2)^(1/2)*sin(S(i))^2*cos(S(i)))/3*(1+cos(S(i))^2)^3;
                vsr(i,1)=ur(i,1)/(1+cos(S(i))^2)^(1/2);
            end
            % reference state 2
            [ Xr_N,Sr_N,Ur_N ] = Ref_State_Mat( s0,uR,Hp,dt );
            
            %==========================================================================%
            % cost function
            %==========================================================================%
            cost = 0;
            
            for i=1:1:Hp
                
                Xaug_r1=[xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1);S(i)];
                Xaug=[X(:,i);S(i)];
                cost=cost+alpha*(Xaug_r1-Xaug)'*Q1*(Xaug_r1-Xaug);
                
                Xaug_r2=[Xr_N(:,i);Sr_N(i)];
                cost=cost+(1-alpha)*(Xaug_r2-Xaug)'*Q2*(Xaug_r2-Xaug);
            end
            
            for i=1:1:Hu
                Ur1=[Fur(i,1);Fvr(i,1);Frr(i,1);vsr(i,1)];
                cost=cost+alpha*(Ur1-U(:,i))'*R1*(Ur1-U(:,i));
                cost=cost+(1-alpha)*(Ur_N(:,i)-U(:,i))'*R2*(Ur_N(:,i)-U(:,i));
            end
            
            Xaugf_r1=[xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1);S(Hp)];
            Xaugf=[X(:,Hp);S(Hp)];
            cost=cost+alpha*(Xaugf_r1-Xaugf)'*Qf1*(Xaugf_r1-Xaugf);
            
            Xaugf_r2=[Xr_N(:,Hp);Sr_N(Hp)];
            cost=cost+(1-alpha)*(Xaugf_r2-Xaugf)'*Qf2*(Xaugf_r2-Xaugf);
        end
        
        %         function [y,ceq] = mpc_contraints(obj)
        %             ceq=[];
        %             y = [];
        %         end
        
        
        
        function u = calc_control(obj,uR,alpha,s0,X0,u0,dt)
            options = optimset('Algorithm','sqp');
            u = fmincon(@(u) obj.mpc_cost(uR,alpha,u,X0,s0,dt),u0,[],[],[],[],obj.lowerbound,obj.upperbound,[],options);
        end
        
        
        %%%%%%%%%%%%%%%%%%
        
        function get_params(obj)
            fprintf('weights: \n')
            disp(obj.weights)
        end
    end
end


