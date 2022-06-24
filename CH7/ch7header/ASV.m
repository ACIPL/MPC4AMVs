classdef ASV < handle 
    properties
        Coef 
        Ndof = 1
        DeducedCoef
        X
        U
    end
    methods
        function obj = ASV(coef,ndof,X0,U0)
            if nargin ==  4
                obj.Coef = coef;
                obj.Ndof = ndof;
                obj = calc_deduced_coef(obj);
                obj.X = X0;
                obj.U = U0;
            else
                errID = 'myComponent:inputError';
                msg = 'Four input arguments (coef, ndof, X0, U0) should be provided.';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        
        function obj = calc_deduced_coef(obj)
            % Currently we only use 1DOF ASV model, but the 3 DOF model
            % should be conveniently accommodated
            if obj.Ndof == 1
                assert(length(obj.Coef) == 4,'Number of model coeffients is incorrect: For 3 DOF AUV model we need 11 parameters.')
                m = obj.Coef(1);
                X_udot = obj.Coef(2);
                Mx = m-X_udot;
                obj.DeducedCoef = [Mx];
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        
        function Xplus = dynamics_discrete(obj,X,U,dt)
            Mx = obj.DeducedCoef(1);
            Xu = obj.Coef(3);
            Du = obj.Coef(4);

            x = X(1);
            u = X(2);

            x_dot=u;
            u_dot=-(Xu/Mx)*u-(Du/Mx)*u*abs(u)+U/Mx;
            X_dot=[x_dot;u_dot];

            Xplus=X+X_dot*dt;
        end
        
        function obj = advance(obj,U,W,dt)
            if obj.Ndof == 1
                disturbed_U = U + W;
                obj.X = obj.dynamics_discrete(obj.X,disturbed_U,dt);
                obj.U = U;
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        
        function Xe  = ErrorState(obj, X, Xr)
            if obj.Ndof == 1
                x = X(1);
                u = X(2);
                xr = Xr(1);
                ur = Xr(2);
                xe = xr-x;
                ue = u-ur;
                
                Xe = [xe;ue];
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end
            
        end
        function f1 = F(obj, X, Xe, P )
            % P = [xR;uR;uRdot];
            % System coefficients
            %==============================================================
            Mx = obj.DeducedCoef(1);
            Xu = obj.Coef(3);
            Du = obj.Coef(4);
            %==============================================================
            x = X(1);
            u = X(2);
            
            xe = Xe(1);
            ue = Xe(2);
            
            xr = P(1);
            ur = P(2);
            ur_dot = P(3);
            
            %==============================================================
            f1 = Xu*u+Du*u*abs(u)+Mx*(ur_dot+ur);
        end
    end
end
