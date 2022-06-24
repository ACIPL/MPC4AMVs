classdef AUV < handle %Model
    properties
        Coef = [116.0; 13.1;-167.6;-477.2;-15.9;26.9;35.8;3.5;241.3;503.8;76.9]
        Ndof = 3;
        DeducedCoef = [283.6; 593.2; 29.0];
        X
        U
    end
    methods
        function obj = AUV(coef,ndof,X0,U0)
            if nargin ==  4
                % Currently 3DOF model is supported and should be extended
                % in the future
%                 if ndof == 3
%                     assert(length(X0) == 6,'For 3 DOF AUV model X = [x;y;psi;u;v;r].');
%                     assert(length(coef) == 11,'coef = [m; Iz; X_udot; Y_vdot; N_rdot; Xu; Yv; Nr; Du; Dv; Dr]')
%                 end
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
                Mx = m-X_udot;
                My = m-Y_vdot;
                Mpsi = Iz-N_rdot;
                obj.DeducedCoef = [Mx;My;Mpsi];
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        
        function Xplus = dynamics_discrete(obj,X,U,dt)
            Mx = obj.DeducedCoef(1);
            My = obj.DeducedCoef(2);
            Mpsi = obj.DeducedCoef(3);
            Xu = obj.Coef(6);
            Yv = obj.Coef(7);
            Nr = obj.Coef(8);
            Du = obj.Coef(9);
            Dv = obj.Coef(10);
            Dr = obj.Coef(11);
%             Du=100.3*(1+err_model);

            x=X(1);
            y=X(2);
            psi=X(3);
            u=X(4);
            v=X(5);
            r=X(6);

            Fu=U(1);
            Fv=U(2);
            Fr=U(3);

            x_dot=u*cos(psi)-v*sin(psi);
            y_dot=u*sin(psi)+v*cos(psi);
            psi_dot=r;
            u_dot=(My/Mx)*v*r-(Xu/Mx)*u-(Du/Mx)*u*abs(u)+Fu/Mx;
            v_dot=-(Mx/My)*u*r-(Yv/My)*v-(Dv/My)*v*abs(v)+Fv/My;
            r_dot=((Mx-My)/Mpsi)*u*v-(Nr/Mpsi)*r-(Dr/Mpsi)*r*abs(r)+Fr/Mpsi;

            X_dot=[x_dot;y_dot;psi_dot;u_dot;v_dot;r_dot];

            Xplus=X+X_dot*dt;
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

        function  hatXplus = dynamics_discrete_ESO(obj,X,U,hatX,K,dt)
            Mx = obj.DeducedCoef(1);
            My = obj.DeducedCoef(2);
            Mpsi = obj.DeducedCoef(3);
            Xu = obj.Coef(6);
            Yv = obj.Coef(7);
            Nr = obj.Coef(8);
            Du = obj.Coef(9);
            Dv = obj.Coef(10);
            Dr = obj.Coef(11);
            M = diag([Mx,My,Mpsi]);
            u=X(4);
            v=X(5);
            r=X(6);
            % For the ESO design
            e1 = X(1:3) - hatX(1:3);
            e2 = X(4:6) - hatX(4:6);
            K1 = K{1};  K2 = K{2}; K3 = K{3};
            hatXplus  = zeros(9,1);
            hatXplus(1) = hatX(4)*cos(hatX(3)) - hatX(5)*sin(hatX(3))+K1*e1(1);
            hatXplus(2) = hatX(4)*sin(hatX(3)) + hatX(5)*cos(hatX(3))+K1*e1(2);
            hatXplus(3) = hatX(6)+K1*e1(3);
            
            hatXplus(4) = (My/Mx)*v*r-(Xu/Mx)*u-(Du/Mx)*u*abs(u)+U(1)/Mx+hatX(7);
            hatXplus(5) = -(Mx/My)*u*r-(Yv/My)*v-(Dv/My)*v*abs(v)+U(2)/My+hatX(8);
            hatXplus(6) = ((Mx-My)/Mpsi)*u*v-(Nr/Mpsi)*r-(Dr/Mpsi)*r*abs(r)+U(3)/Mpsi+hatX(9);
            hatXplus(4:6) = hatXplus(4:6) +  K2 * e1;
            hatXplus(7:9) = K3*e1;
            hatXplus = hatX + hatXplus * dt;
        end
        
        function obj = advance_ESO(obj,U,W,dt)
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
        
        function Xe  = ErrorState(obj, X, Xr)
            if obj.Ndof == 3
                x = X(1);
                y = X(2);
                psi = X(3);
                u = X(4);
                v = X(5);
                r = X(6);
                
                xr = Xr(1);
                yr = Xr(2);
                psi_r = Xr(3);
                ur = Xr(4);
                vr = Xr(5);
                rr = Xr(6);
                
                xe = (xr-x)*cos(psi)+(yr-y)*sin(psi);
                ye = -(xr-x)*sin(psi)+(yr-y)*cos(psi);
                psi_e = psi_r-psi;
                ue = u-ur*cos(psi_e)+vr*sin(psi_e);
                ve = v-ur*sin(psi_e)-vr*cos(psi_e);
                re = r-rr;
                
                Xe = [xe;ye;psi_e;ue;ve;re];
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end
            
        end
        function [ f1, f2, f3 ] = F(obj, X,Xe,P )
            % P = [xR;yR;psiR;uR;vR;rR;uRdot;vRdot;rRdot];
            % System coefficients
            %==========================================================================
            Mx = obj.DeducedCoef(1);
            My = obj.DeducedCoef(2);
            Mpsi = obj.DeducedCoef(3);
            Xu = obj.Coef(6);
            Yv = obj.Coef(7);
            Nr = obj.Coef(8);
            Du = obj.Coef(9);
            Dv = obj.Coef(10);
            Dr = obj.Coef(11);
            %==========================================================================
            x = X(1);
            y = X(2);
            psi = X(3);
            u = X(4);
            v = X(5);
            r = X(6);
            
            xe = Xe(1);
            ye = Xe(2);
            psi_e = Xe(3);
            ue = Xe(4);
            ve = Xe(5);
            re = Xe(6);
            
            xr = P(1);
            yr = P(2);
            psi_r = P(3);
            ur = P(4);
            vr = P(5);
            rr = P(6);
            ur_dot = P(7);
            vr_dot = P(8);
            rr_dot = P(9);
            
            %==========================================================================
            
            f1 = -My*v*r+Xu*u+Du*u*abs(u)+Mx*(ur_dot*cos(psi_e)+ur*sin(psi_e)*re-vr_dot*sin(psi_e)+vr*cos(psi_e)*re);
            f2 = Mx*u*r+Yv*v+Dv*v*abs(v)+My*(ur_dot*sin(psi_e)-ur*cos(psi_e)*re+vr_dot*cos(psi_e)+vr*sin(psi_e)*re);
            f3 = (My-Mx)*u*v+Nr*r+Dr*r*abs(r)+Mpsi*rr_dot;
        end
        
        
    end
end
