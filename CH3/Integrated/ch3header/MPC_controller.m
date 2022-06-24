classdef MPC_controller < handle
    properties
        N
        model %Model = AUV(ones(11,1),3,zeros(6,1),zeros(3,1))
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
        
        function cost = mpc_cost(obj,u,X0,P,dt)
            Mx = obj.model.DeducedCoef(1);
            My = obj.model.DeducedCoef(2);
            Mpsi = obj.model.DeducedCoef(3);
            
            Hp = obj.N;
            Q = obj.weights{1};
            R = obj.weights{2};
            Qf = obj.weights{3};
            
            nu = length(u);
            nu = nu/Hp;
            nx = length(X0);
            
            U = zeros(nu,Hp);
            X = zeros(nx,Hp);
            
            xr=zeros(Hp,1);
            yr=zeros(Hp,1);
            psi_r=zeros(Hp,1);
            ur=zeros(Hp,1);
            vr=zeros(Hp,1);
            rr=zeros(Hp,1);
            urdot=zeros(Hp,1);
            vrdot=zeros(Hp,1);
            rrdot=zeros(Hp,1);
            
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
%             Xplus = AUV_Dynamics_discrete( X0,U(:,1),dt );
            X(:,1) = Xplus;
            
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,U(:,i),dt );
%                 Xplus = AUV_Dynamics_discrete(Xplus,U(:,i),dt);
                X(:,i)= Xplus;
            end
            
            %==========================================================================%
            % reference states
            %==========================================================================%
            for i=1:1:Hp
                
                xr(i,1)=P(1,i);
                yr(i,1)=P(2,i);
                psi_r(i,1)=P(3,i);
                ur(i,1)=P(4,i);
                vr(i,1)=P(5,i);
                rr(i,1)=P(6,i);
                
                urdot(i,1)=P(7,i);
                vrdot(i,1)=P(8,i);
                rrdot(i,1)=P(9,i);
                
            end
            
            
            %==========================================================================%
            % cost function
            %==========================================================================%
            cost=0;
            
            for i=1:1:Hp
                
                Xr=[xr(i,1);yr(i,1);psi_r(i,1);ur(i,1);vr(i,1);rr(i,1)];
                Xe = obj.model.ErrorState( X(:,i), Xr );
                [f1, f2, f3] = obj.model.F( X(:,i),Xe,P(:,i) );
                tau_u = (U(1,i)-f1)/Mx;
                tau_v = (U(2,i)-f2)/My;
                tau_r = (U(3,i)-f3)/Mpsi;
                Te = [tau_u;tau_v;tau_r];
                cost=cost+Xe'*Q*Xe+Te'*R*Te;
                
            end
            
            Xr=[xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
            Xe = obj.model.ErrorState( X(:,Hp), Xr );
            cost=cost+Xe'*Qf*Xe;
        end
        
        function [ y,ceq ] = mpc_contraints( obj,u,X0,P,dt)
            
            ceq=[];
            
            %==========================================================================
            % MPC
            Hp = obj.N;
            nx=length(X0);
            nu=length(u); %[U0,U1,..U_N-1]
            nu=nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
            U=zeros(nu,Hp);
            
            xr=zeros(Hp,1);
            yr=zeros(Hp,1);
            psi_r=zeros(Hp,1);
            ur=zeros(Hp,1);
            vr=zeros(Hp,1);
            rr=zeros(Hp,1);
            urdot=zeros(Hp,1);
            vrdot=zeros(Hp,1);
            rrdot=zeros(Hp,1);
            
            
            %==============================================================
            % partition of u | Ui=[Fu_i, Fv_i, Fr_i]'
            
            for i=1:1:Hp
                for j=1:1:nu
                    
                    U(j,i)=u((i-1)*nu+j,1);
                    
                end
            end
            
            %==============================================================
            % reference state
            
            for i=1:1:Hp
                
                xr(i,1)=P(1,i);
                yr(i,1)=P(2,i);
                psi_r(i,1)=P(3,i);
                ur(i,1)=P(4,i);
                vr(i,1)=P(5,i);
                rr(i,1)=P(6,i);
                
                urdot(i,1)=P(7,i);
                vrdot(i,1)=P(8,i);
                rrdot(i,1)=P(9,i);
                
            end
            
            %==========================================================================
            %state prediction
            
            Xall=zeros(nx,Hp+1);
            Xplus=X0;
            Xall(:,1)=Xplus;
            
            for i=1:1:Hp-1
%                 Xplus = AUV_Dynamics_discrete(Xplus,U(:,i),dt);
                Xplus = obj.model.dynamics_discrete(Xplus,U(:,i),dt);
                Xall(:,i+1)=Xplus;
                
            end
            
            XN = Xall(:,Hp);
            XrN = [ xr(Hp,1);yr(Hp,1);psi_r(Hp,1);ur(Hp,1);vr(Hp,1);rr(Hp,1)];
            XeN = obj.model.ErrorState( XN, XrN );
            xeN = XeN(1);
            yeN = XeN(2);
            psi_eN = XeN(3);
            ueN = XeN(4);
            veN = XeN(5);
            reN = XeN(6);
            
            y1 = abs(xeN) - abs(ueN);
            y2 = abs(yeN) - abs(veN);
            y3 = abs(psi_eN) - abs(reN);
            y4 = -xeN*ueN;
            y5 = -yeN*veN;
            y6 = -psi_eN*reN;
            
            y=[y1;y2;y3;y4;y5;y6];
        end
        
        function u = calc_control(obj,u0,X0,P,dt)
             options = optimset('Algorithm','sqp');
            u = fmincon(@(u) obj.mpc_cost(u,X0,P,dt),u0,[],[],...
                [],[],obj.lowerbound,obj.upperbound,...
                @(u)obj.mpc_contraints(u,X0,P,dt),options);
        end

        %%%%%%%%%%%%%%%%%%
        
        function get_params(obj)
            fprintf('weights: \n')
            disp(obj.weights)
        end
    end
end


