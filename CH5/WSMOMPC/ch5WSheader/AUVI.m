classdef AUVI < Model
    properties
        Coef = [116.0; 13.1;-167.6;-477.2;-15.9;26.9;35.8;3.5;241.3;503.8;76.9]
        Ndof = 3;
        DeducedCoef = [283.6; 593.2; 29.0]
    end
    methods
        function obj = AUVI(coef,ndof,X0,U0)
            if nargin == 4
                % Currently 3DOF model is supported and should be extended
                % in the future
                if ndof == 3
                    assert(length(X0) == 6,'For 3 DOF AUV model X = [x;y;psi;u;v;r].');
                    assert(length(coef) == 11,'coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr]')     
                end       
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
            % Currently we only use 3DOF AUV model, but the 6 DOF model
            % should be conveniently accommodated
            if obj.Ndof == 3 
                assert(length(obj.Coef) == 11,'Number of model coeffients is incorrect: For 3 DOF AUV model we need 11 parameters.')        
                m = obj.Coef(1);
                Iz = obj.Coef(2);
                X_udot = obj.Coef(3);
                Y_vdot = obj.Coef(4);
                N_rdot = obj.Coef(5);
                Mudot = m-X_udot;
                Mvdot = m-Y_vdot;
                Mrdot = Iz-N_rdot;
                obj.DeducedCoef = [Mudot;Mvdot;Mrdot];
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end   
        end
        
        function Xplus = dynamics_discrete(obj,X,U,dt)
                Mudot = obj.DeducedCoef(1);
                Mvdot = obj.DeducedCoef(2);
                Mrdot = obj.DeducedCoef(3);
                Xu = obj.Coef(6);
                Yv = obj.Coef(7);
                Nr = obj.Coef(8);
                Du = obj.Coef(9);
                Dv = obj.Coef(10);
                Dr = obj.Coef(11);
                tau = U;
                Fu = tau(1);
                Fv = tau(2);
                Fr = tau(3);
                psi = X(3); 
                u = X(4); 
                v = X(5); 
                r = X(6); 
                x_dot = u*cos(psi)-v*sin(psi);
                y_dot = u*sin(psi)+v*cos(psi);
                psi_dot = r;
                u_dot = (Mvdot/Mudot)*v*r-(Xu/Mudot)*u-(Du/Mudot)*u*abs(u)+ Fu/Mudot;
                v_dot = -(Mudot/Mvdot)*u*r-(Yv/Mvdot)*v-(Dv/Mvdot)*v*abs(v)+ Fv/Mvdot;
                r_dot = ((Mudot-Mvdot)/Mrdot)*u*v-(Nr/Mrdot)*r-(Dr/Mrdot)*r*abs(r)+ Fr/Mrdot;
                X_dot = [x_dot;y_dot;psi_dot;u_dot;v_dot;r_dot];
                Xplus = X + X_dot*dt; 
        end
        
        function obj = advance(obj,U,W,dt)
            if obj.Ndof == 3 
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

    end
end
