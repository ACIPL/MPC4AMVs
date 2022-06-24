classdef LMPC_controller < Controller
    properties
        N
        model Model = AUV(ones(11,1),3,eye(3),zeros(6,1),zeros(3,1))
        auxiliary_controller Auxuliary_Controller = nonlinearPD_controller(eye(3),eye(3))
        weights
        upperbound
        lowerbound       
    end
    
    methods
        function obj = LMPC_controller(N,model,auxiliary_controller,weights,upperbound,lowerbound)
            obj.N = N;
            obj.model = model;
            obj.auxiliary_controller = auxiliary_controller;
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
        
        function cost = lmpc_cost(obj,u,X0,P,dt)
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

            T1 = zeros(Hp,1);
            T2 = zeros(Hp,1);
            T3 = zeros(Hp,1);
            T4 = zeros(Hp,1);

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
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,U(:,i),dt );
                X(:,i)= Xplus;
            end

            %==========================================================================%
            % reference state
            %==========================================================================%
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
                T1(i,1) = P(10,i);
                T2(i,1) = P(11,i);
                T3(i,1) = P(12,i);
                T4(i,1) = P(13,i);
            end    

            %==========================================================================%
            % cost function
            %==========================================================================%
            cost = 0;

            for i=1:1:Hp
                Xr = [xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1)];
                cost = cost + (Xr-X(:,i))'*Q*(Xr-X(:,i));
            end

            for i=1:1:Hu
                cost = cost + U(:,i)'*R*U(:,i);
            end

            Xr = [xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
            cost = cost + (Xr-X(:,Hp))'*Qf*(Xr-X(:,Hp));
        end
        
        function [y,ceq] = lmpc_contraints(obj,u,X0,B,P)
            ceq=[];
            Hp = obj.N;
            nu = length(u); %[U0,U1,..U_N-1]
            nu = nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
            U = zeros(nu,Hp); 
            
            %==========================================================================
            % partition of u | Ui=[Fu_i, Fv_i, Fr_i]'
            %==========================================================================
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i)=u((i-1)*nu+j,1);
                end
            end  
            Tau = B*U(:,1);    
            eta_d = P(1:3);
            Tau_auxiliary = obj.auxiliary_controller.calc_control(P,X0); 
            Vdot_auxiliary = obj.auxiliary_controller.lyapunov_derivative(X0,Tau_auxiliary,eta_d); 
            y = obj.auxiliary_controller.lyapunov_derivative(X0,Tau,eta_d) - Vdot_auxiliary;     
%             y =[];
        end
        
     

        function u = calc_control(obj,P,X0,B,u0,dt)
            options = optimset('Algorithm','sqp');
            u = fmincon(@(u) obj.lmpc_cost(u,X0,P,dt),u0,[],[],[],[],obj.lowerbound,obj.upperbound,@(u) obj.lmpc_contraints( u,X0,B,P),options);
        end
       
        function u0 = calc_initial_guess(obj,P,X0,CtrlMat_pinv,dt)
            Hp = obj.N;
            [row,col] = size(P);
            assert(Hp == col,'Reference is incompatible with the prediction horizon.')
            nu = length(obj.model.U);
            u0 = zeros(nu*Hp,1);
            X0_ = X0;
            for j = 1:1:Hp
                U_auxiliary = obj.auxiliary_controller.calc_control(P(:,j),X0_);
                U_auxiliary_ = CtrlMat_pinv * U_auxiliary; 
                u0(nu*(j-1)+1:nu*j) = U_auxiliary_;
                X0_ = obj.model.dynamics_discrete(X0_,U_auxiliary_,dt);
            end            
        end
               
        %%%%%%%%%%%%%%%%%% 
        
        function get_params(obj)
            fprintf('weights: \n')
            disp(obj.weights)     
        end
    end
end


