classdef RDMPC_controller < Controller
    properties
        N
        model %Model = AUV(ones(11,1),3,zeros(6,1),zeros(3,1))
        weights
        upperbound
        lowerbound
    end
    
    methods
        function obj = RDMPC_controller(N,model,weights,upperbound,lowerbound)
            obj.N = N;
            obj.model = model;
            obj.weights = weights;
            obj.calc_upperbound(upperbound);
            obj.calc_lowerbound(lowerbound);
        end
        
        function obj = calc_upperbound(obj,upperbound)
            if ~isempty(upperbound)
                m = length(upperbound);
                U_max = zeros(m*obj.N,1);
                for i=1:1:obj.N
                    U_max(1+m*(i-1):m*i,1) = upperbound;
                end
                obj.upperbound = U_max;
            else
                obj.upperbound = [];
            end
        end
        
        function obj = calc_lowerbound(obj,lowerbound)
            if ~isempty(lowerbound)
                m = length(lowerbound);
                U_min = zeros(m*obj.N,1);
                for i=1:1:obj.N
                    U_min(1+m*(i-1):m*i,1) = lowerbound;
                end
                obj.lowerbound = U_min;
            else
                obj.lowerbound = [];
            end
        end
        
        function cost = rdmpc_cost(obj,u,X0,P,X_neighbors,dt,dr)
            Hp = obj.N;
            Q = obj.weights{1};
            R = obj.weights{2};
            Qf = obj.weights{3};
            F = obj.weights{4};
            H = obj.weights{5};
            nu = length(u);
            nu = nu/Hp;
            nx = length(X0);
            Hu = Hp;
            U = zeros(nu,Hp);
            X = zeros(nx,Hp);
            
            xr = zeros(Hp,1);
            yr = zeros(Hp,1);
            psi_r = zeros(Hp,1);
            ur = zeros(Hp,1);
            vr = zeros(Hp,1);
            rr = zeros(Hp,1);
            Fur = zeros(Hp,1);
            Fvr = zeros(Hp,1);
            Frr = zeros(Hp,1);
            
            % partition of u
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i) = u((i-1)*nu+j,1);
                end
            end
            
            % state prediction
            Xplus = obj.model.dynamics_discrete( X0,U(:,1),dt );
            X(:,1) = Xplus;
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,U(:,i),dt );
                X(:,i)= Xplus;
            end
            
            % reference state
            for i=1:1:Hp
                xr(i) = P(1,i);
                ur(i) = P(2,i);
            end
            
            % cost function
            cost = 0;
            if length(X_neighbors) == 2
                for i=1:1:Hp - 1
                    Xr = [xr(i);ur(i)];
                    cost = cost + (Xr-X(:,i))'*Q*(Xr-X(:,i))...
                        + (X(:,i)-dr-X_neighbors{1}(:,i))'*F*(X(:,i)-dr-X_neighbors{1}(:,i))...
                        + (X(:,i)-dr-X_neighbors{2}(:,i))'*H*(X(:,i)-dr-X_neighbors{2}(:,i));
                end
            else
                for i=1:1:Hp - 1
                    Xr = [xr(i);ur(i)];
                    cost = cost + (Xr-X(:,i))'*Q*(Xr-X(:,i)) ...
                        + (X(:,i) - dr -X_neighbors{1}(:,i))'*F*(X(:,i) - dr -X_neighbors{1}(:,i));
                end
            end
            
            for i=1:1:Hu
                cost = cost + U(i)'*R*U(i);
            end
            
            cost = cost + (P(:,Hp)-X(:,Hp))'*Qf*(P(:,Hp)-X(:,Hp));
        end
        
        
        function [c,ceq] = rdmpc_contraints(obj,u,X0, X_neighbors,P,dt,dr)
            ceq=[];
            c = [];
            Hp = obj.N;
            nu = length(u); %[U0,U1,..U_N-1]
            nu = nu/Hp; % Ui=[Fu_i]'
            nx = length(X0);
            U = zeros(nu,Hp);
            X = zeros(nx,Hp);
            
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i)=u((i-1)*nu+j,1);
                end
            end
            nei_num = length(X_neighbors);
            
            Xplus = obj.model.dynamics_discrete( X0,U(:,1),dt );
            X(:,1) = Xplus;
            if nei_num == 1
                for i=1:1:Hp
                    if i < Hp
                        Xplus = obj.model.dynamics_discrete( Xplus,U(:,i),dt );
                        X(:,i)= Xplus;
                        y = norm(X(:,i) -dr - X_neighbors{1}(:,i)) - 0.5;
                        c = [c y];
                    elseif i == Hp
                        y = norm(X(:,Hp)- P(1:2,Hp)) - 0.2;
                        c = [c y];
                    end
                end
            elseif nei_num == 2
                for i=1:1:Hp
                    if i < Hp
                        Xplus = obj.model.dynamics_discrete( Xplus,U(:,i),dt );
                        X(:,i)= Xplus;
                        y(1) = norm(X(:,i) - dr - X_neighbors{1}(:,i)) - 0.5;
                        y(2) = X(1,i) - dr(1) -2.5 - X_neighbors{2}(1,i) + 0.5;
                        c = [c y];
                    elseif i == Hp
                        y = norm(X(:,Hp)- P(1:2,Hp)) - 0.2;
                        c = [c y];
                    end
                end
            end
        end
        
        % general control
        function [u, X] = calc_control(obj,P, X_neighbors, X0, u0, dt, dr)
            options = optimset('Algorithm','sqp');
            % Calculate the MPC control input
            u = fmincon(@(u) obj.rdmpc_cost(u,X0,P,X_neighbors,dt,dr),u0,[],[],[],[],...
                obj.lowerbound, obj.upperbound, ...
                @(u) obj.rdmpc_contraints( u,X0,X_neighbors,P,dt,dr),options);
            
            % Construct the broadcast information
            Hp = obj.N;
            Xplus = obj.model.dynamics_discrete( X0, u(1),dt );
            X(:,1) = Xplus - dr(1:2);
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,u(i),dt );
                X(:,i)= Xplus - dr(1:2);
            end
        end
        
        function u0 = calc_initial_guess(obj,P,X0,dt)
            Hp = obj.N;
            nu = length(obj.model.U);
            u0 = zeros(nu*Hp,1);
            X0_ = X0;
            for j = 1:1:Hp
                U_auxiliary = obj.auxiliary_controller.calc_control(P(:,j),X0_);
                u0(nu*(j-1)+1:nu*j) = U_auxiliary;
                X0_ = obj.model.dynamics_discrete(X0_,U_auxiliary,dt);
            end
        end
        
        function get_params(obj)
            fprintf('weights: \n')
            disp(obj.weights)
        end
    end
end


