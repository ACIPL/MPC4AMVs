classdef DLMPC_controller < Controller
    properties
        N
        model %Model = AUV(ones(11,1),3,zeros(6,1),zeros(3,1))
        auxiliary_controller Auxuliary_Controller = nonlinearBS_controller(eye(3),eye(3))
        weights
        upperbound
        lowerbound
    end
    
    methods
        function obj = DLMPC_controller(N,model,controller,weights,upperbound,lowerbound)
            obj.N = N;
            obj.model = model;
            obj.auxiliary_controller = controller;
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
        
        function cost = dlmpc_cost(obj,u,X0,P,X_neighbors,dt,dr)
            Hp = obj.N;
            Q = obj.weights{1};
            R = obj.weights{2};
            Qf = obj.weights{3};
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
                xr(i,1) = P(1,i);
                yr(i,1) = P(2,i);
                psi_r(i,1) = P(3,i);
                ur(i,1) = P(4,i);
                vr(i,1) = P(5,i);
                rr(i,1) = P(6,i);
                Fur(i,1) = P(7,i);
                Fvr(i,1) = P(8,i);
                Frr(i,1) = P(9,i);
            end
            
            % cost function
            cost = 0;
            for i=1:1:Hp - 1
                Xr = [xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1)];
                cost = cost + (Xr-X(:,i))'*Q*(Xr-X(:,i)) ...
                    + (X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))...
                    + (X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i));
            end
            
            for i=1:1:Hu
                cost = cost + U(:,i)'*R*U(:,i);
            end
            
            cost = cost + (Xr-X(:,Hp))'*Qf*(Xr-X(:,Hp)) ...
                + (X(1:2,Hp)-dr(1:2)-X_neighbors{1}(1:2,Hp))'*Qf(1:2,1:2)*(X(1:2,Hp)...
                -dr(1:2)-X_neighbors{1}(1:2,Hp))...
                + (X(1:2,Hp)-dr(1:2)-X_neighbors{2}(1:2,Hp))'*Qf(1:2,1:2)*(X(1:2,Hp)...
                -dr(1:2)-X_neighbors{2}(1:2,Hp));
        end
        
        function cost = dlmpc_cost_avoidance(obj,u,X0,P,X_neighbors,dt,dr)
            Hp = obj.N;
            Q = obj.weights{1};
            R = obj.weights{2};
            Qf = obj.weights{3};
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
                xr(i,1) = P(1,i);
                yr(i,1) = P(2,i);
                psi_r(i,1) = P(3,i);
                ur(i,1) = P(4,i);
                vr(i,1) = P(5,i);
                rr(i,1) = P(6,i);
                Fur(i,1) = P(7,i);
                Fvr(i,1) = P(8,i);
                Frr(i,1) = P(9,i);
            end
            
            % cost function
            cost = 0;
            d1r = [0;1]; d2r = [0;-1]; d3r = [-1;0];
            for i=1:1:Hp - 1
                Xr = [xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1)];
                if dr(1) == d1r(1) && dr(2) == d1r(2)
                    cost = cost + (Xr-X(:,i))'*Q*(Xr-X(:,i)) ...
                        + (X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))...
                        + (X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i))...
                        + 1e4/(1+exp(10*norm(X(1:2,i) - X_neighbors{1}(1:2,i) - d2r) - 6))...
                        + 1e4/(1+exp(10*norm(X(1:2,i) - X_neighbors{2}(1:2,i) - d3r) - 6));
                elseif dr(1) == d2r(1) && dr(2) == d2r(2)
                    cost = cost + (Xr-X(:,i))'*Q*(Xr-X(:,i)) ...
                        + (X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))...
                        + (X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i))...
                        + 1e4/(1+exp(10*norm(X(1:2,i) - X_neighbors{1}(1:2,i) - d1r) - 6))...
                        + 1e4/(1+exp(10*norm(X(1:2,i) - X_neighbors{2}(1:2,i) - d2r) - 6));
                elseif dr(1) == d3r(1) && dr(2) == d3r(2)
                    cost = cost + (Xr-X(:,i))'*Q*(Xr-X(:,i)) ...
                        + (X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{1}(1:2,i))...
                        + (X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i))'*Q(1:2,1:2)*(X(1:2,i)-dr(1:2)-X_neighbors{2}(1:2,i))...
                        + 1e4/(1+exp(10*norm(X(1:2,i) - X_neighbors{1}(1:2,i) - d1r) - 6))...
                        + 1e4/(1+exp(10*norm(X(1:2,i) - X_neighbors{2}(1:2,i) - d2r) - 6));
                end
            end
            
            for i=1:1:Hu
                cost = cost + U(:,i)'*R*U(:,i);
            end
            
            cost = cost + (Xr-X(:,Hp))'*Qf*(Xr-X(:,Hp)) ...
                + (X(:,Hp)-dr(1:6)-X_neighbors{1}(:,Hp))'*Qf*(X(:,Hp)...
                -dr(1:6)-X_neighbors{1}(:,Hp))...
                + (X(:,Hp)-dr(1:6)-X_neighbors{2}(:,Hp))'*Qf*(X(:,Hp)...
                -dr(1:6)-X_neighbors{2}(:,Hp));
        end
        
        function [y,ceq] = dlmpc_contraints(obj,u,X0, P,traj,t,dr)
            ceq=[];
            Hp = obj.N;
            nu = length(u); %[U0,U1,..U_N-1]
            nu = nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
            U = zeros(nu,Hp);
            
            M = diag([obj.model.DeducedCoef]);
            PI=[0,-1,0;1,0,0;0,0,0];
            Lambda = diag([1;1;1]);
            psi=X0(3);
            Kp = obj.auxiliary_controller.Kp;
            Kd = obj.auxiliary_controller.Kd;
            eta=X0(1:3,1);
            vel=X0(4:6,1);
            C=CV(vel);
            D=DV(vel);
            
            Rpsi=Rot(psi);
            Dstar=Rpsi*D*Rpsi';
            
            xR = traj.Xr(t);
            xRdot = traj.Xrdot(t);
            xRddot = traj.Xrddot(t);
            yR = traj.Yr(t);
            yRdot = traj.Yrdot(t);
            yRddot = traj.Yrddot(t);
            psiR = traj.PSIr(t);
            psiRdot = traj.PSIrdot(t);
            psiRddot = traj.PSIrddot(t);
            
            eta_d=[xR;yR;psiR]+dr(1:3);
            eta_d_dot=[xRdot;yRdot;psiRdot];
            eta_d_ddot=[xRddot;yRddot;psiRddot];
            
            tilde_eta=eta-eta_d;
            eta_r_dot=eta_d_dot-Lambda*tilde_eta;
            
            v_r=Rpsi'*eta_r_dot;
            
            eta_dot=Rpsi*vel;
            S=eta_dot-eta_r_dot;
            
            psi_dot=eta_dot(3);
            v_r_dot=-psi_dot*Rpsi'*PI*eta_r_dot+Rpsi'*(eta_d_ddot-...
                Lambda*(eta_dot-eta_d_dot));
            tau=M*v_r_dot+C*v_r+D*v_r-Rpsi'*Kp*tilde_eta-Rpsi'*Kd*S;
            
            V2dot_backstepping=S'*Rpsi*(tau-M*v_r_dot-C*v_r-D*v_r+Rpsi'*...
                Kp*tilde_eta)-S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta;
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i)=u((i-1)*nu+j,1);
                end
            end
            Tau=U(:,1);
            y=S'*Rpsi*(Tau-M*v_r_dot-C*v_r-D*v_r+Rpsi'*Kp*tilde_eta)...
                -S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta-V2dot_backstepping;
            
        end
        
        function [y,ceq] = dlmpc_contraints_ESO(obj,u,X0, hatW,traj,t,dr)
            ceq=[];
            Hp = obj.N;
            nu = length(u); %[U0,U1,..U_N-1]
            nu = nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
            U = zeros(nu,Hp);
            
            M = diag([obj.model.DeducedCoef]);
            PI=[0,-1,0;1,0,0;0,0,0];
            Lambda = diag([1;1;1]);
            psi=X0(3);
            Kp = obj.auxiliary_controller.Kp;
            Kd = obj.auxiliary_controller.Kd;
            eta=X0(1:3,1);
            vel=X0(4:6,1);
            C=CV(vel);
            D=DV(vel);
            
            Rpsi=Rot(psi);
            Dstar=Rpsi*D*Rpsi';
            
            xR = traj.Xr(t);
            xRdot = traj.Xrdot(t);
            xRddot = traj.Xrddot(t);
            yR = traj.Yr(t);
            yRdot = traj.Yrdot(t);
            yRddot = traj.Yrddot(t);
            psiR = traj.PSIr(t);
            psiRdot = traj.PSIrdot(t);
            psiRddot = traj.PSIrddot(t);
            
            eta_d=[xR;yR;psiR]+dr(1:3);
            eta_d_dot=[xRdot;yRdot;psiRdot];
            eta_d_ddot=[xRddot;yRddot;psiRddot];
            
            tilde_eta=eta-eta_d;
            eta_r_dot=eta_d_dot-Lambda*tilde_eta;
            
            v_r=Rpsi'*eta_r_dot;
            
            eta_dot=Rpsi*vel;
            S=eta_dot-eta_r_dot;
            
            psi_dot=eta_dot(3);
            v_r_dot=-psi_dot*Rpsi'*PI*eta_r_dot+Rpsi'*(eta_d_ddot-...
                Lambda*(eta_dot-eta_d_dot));
            tau=M*v_r_dot+C*v_r+D*v_r-Rpsi'*Kp*tilde_eta-Rpsi'*Kd*S -hatW;
            
            V2dot_backstepping=S'*Rpsi*(tau-M*v_r_dot-C*v_r-D*v_r+Rpsi'*...
                Kp*tilde_eta)-S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta;
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i)=u((i-1)*nu+j,1);
                end
            end
            Tau=U(:,1);
            y=S'*Rpsi*(Tau-M*v_r_dot-C*v_r-D*v_r+Rpsi'*Kp*tilde_eta)...
                -S'*Dstar*S-tilde_eta'*Kp*Lambda*tilde_eta-V2dot_backstepping;
            
        end
        % general control
        function [u, X] = calc_control(obj,P,traj, X_neighbors, X0, u0, dt, t, dr)
            options = optimset('Algorithm','sqp');
            % Calculate the MPC control input
            u = fmincon(@(u) obj.dlmpc_cost(u,X0,P,X_neighbors,dt,dr),u0,[],[],[],[],...
                obj.lowerbound, obj.upperbound, ...
                @(u) obj.dlmpc_contraints( u,X0,P,traj,t,dr),options);
            
            % Construct the broadcast information
            Hp = obj.N;
            Xplus = obj.model.dynamics_discrete( X0, u(1:3),dt );
            X(:,1) = Xplus - dr(1:6);
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,u(3*i-2:3*i),dt );
                X(:,i)= Xplus - dr(1:6);
            end
        end
        %DLMPC considering collision avoidance
        function [u, X] = calc_control_avoidance(obj,P,traj, X_neighbors, X0, u0, dt, t, dr)
            options = optimset('Algorithm','sqp');
            % Calculate the MPC control input
            u = fmincon(@(u) obj.dlmpc_cost_avoidance(u,X0,P,X_neighbors,dt,dr),u0,[],[],[],[],...
                obj.lowerbound, obj.upperbound, ...
                @(u) obj.dlmpc_contraints( u,X0,P,traj,t,dr),options);
            
            % Construct the broadcast information
            Hp = obj.N;
            Xplus = obj.model.dynamics_discrete( X0, u(1:3),dt );
            X(:,1) = Xplus - dr(1:6);
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,u(3*i-2:3*i),dt );
                X(:,i)= Xplus - dr(1:6);
            end
                end
        % DLMPC considering disturbances with the ESO
        function [u, X] = calc_control_ESO(obj,P,traj, X_neighbors, X0, u0,hatW, dt, t, dr)
            options = optimset('Algorithm','sqp');
            % Calculate the MPC control input
            u = fmincon(@(u) obj.dlmpc_cost_avoidance(u,X0,P,X_neighbors,dt,dr),u0,[],[],[],[],...
                obj.lowerbound, obj.upperbound, ...
                @(u) obj.dlmpc_contraints_ESO( u,X0,hatW,traj,t,dr),options);
            
            % Construct the broadcast information
            Hp = obj.N;
            Xplus = obj.model.dynamics_discrete( X0, u(1:3),dt );
            X(:,1) = Xplus - dr(1:6);
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,u(3*i-2:3*i),dt );
                X(:,i)= Xplus - dr(1:6);
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


