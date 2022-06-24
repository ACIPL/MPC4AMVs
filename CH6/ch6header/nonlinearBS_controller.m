classdef nonlinearBS_controller < Auxuliary_Controller
    properties
        Kp
        Kd
    end
    
    methods
        function obj = nonlinearBS_controller(Kp,Kd)
            obj.Kp = Kp;
            obj.Kd = Kd;
        end
        
        function [U, para] = calc_control(obj,ref,state)
            
            m=116;
            Iz=13.1;
            X_udot=-167.6;
            Y_vdot=-477.2;
            N_rdot=-15.9;
           
            Mx=m-X_udot;
            My=m-Y_vdot;
            Mpsi=Iz-N_rdot;
            
            eta_d=ref(1:3);
            eta_d_dot=ref(4:6);
            eta_d_ddot=ref(7:9);
            eta = state(1:3,1); % [x;y;psi]
            vel = state(4:6,1); % [u;v;r]
            psi = state(3);
            Rpsi = Rot(psi);
            M=diag([Mx;My;Mpsi]);
            Lambda = diag([1;1;1]);
            tilde_eta=eta-eta_d;
            eta_r_dot=eta_d_dot-Lambda*tilde_eta;
            
            v_r=Rpsi'*eta_r_dot;
            
            eta_dot=Rpsi*vel;
            S=eta_dot-eta_r_dot;
            
            psi_dot=eta_dot(3);
            PI=[0,-1,0;1,0,0;0,0,0];
            
            v_r_dot=-psi_dot*Rpsi'*PI*eta_r_dot+Rpsi'*(eta_d_ddot-Lambda*(eta_dot-eta_d_dot));
            
            C=CV(vel);
            D=DV(vel);
            g = 0;
            %   Backstepping
            U=M*v_r_dot+C*v_r+D*v_r+g-Rpsi'*obj.Kp*tilde_eta-Rpsi'*obj.Kd*S;
            para = [v_r_dot, v_r, tilde_eta, S];
        end
        
     
        
        function get_params(obj)
            fprintf('Kp: \n')
            disp(obj.Kp)
            fprintf('Kd: \n')
            disp(obj.Kd)
        end
        
        function Vdot = lyapunov_derivative(obj,X,ref)
%             eta=X(1:3);
            vel=X(4:6);
%             tilde_eta = eta - eta_d;
%             Vdot = vel'*(tau-C*vel-D*vel+Rpsi'*obj.Kp*tilde_eta);
            
            M = diag([obj.model.DeducedCoef]);
%             PI=[0,-1,0;1,0,0;0,0,0];
            Lambda = diag([1;1;1]);
            psi=X(3);
            Kp1 = obj.Kp;
            Kd1 = obj.Kd;

            C=CV(vel);
            D=DV(vel);
            
            Rpsi=Rot(psi);
            Dstar=Rpsi*D*Rpsi';
            
            
%             xR = traj.Xr(t);
%             xRdot = traj.Xrdot(t);
%             xRddot = traj.Xrddot(t);
%             yR = traj.Yr(t);
%             yRdot = traj.Yrdot(t);
%             yRddot = traj.Yrddot(t);
%             psiR = traj.PSIr(t);
%             psiRdot = traj.PSIrdot(t);
%             psiRddot = traj.PSIrddot(t);
%             
%             eta_d=[xR;yR;psiR];
%             eta_d_dot=[xRdot;yRdot;psiRdot];
%             eta_d_ddot=[xRddot;yRddot;psiRddot];
%             
%             tilde_eta=eta-eta_d;
%             eta_r_dot=eta_d_dot-Lambda*tilde_eta;
%             
%             v_r=Rpsi'*eta_r_dot;
%             
%             eta_dot=Rpsi*vel;
%             S=eta_dot-eta_r_dot;
%             
%             psi_dot=eta_dot(3);
%             v_r_dot=-psi_dot*Rpsi'*PI*eta_r_dot+Rpsi'*(eta_d_ddot-Lambda*(eta_dot-eta_d_dot));
%             tau=M*v_r_dot+C*v_r+D*v_r-Rpsi'*Kp1*tilde_eta-Rpsi'*Kd1*S;
            [tau,para] = obj.calc_control(obj,ref,X);
            v_r_dot = para(1); 
            v_r = para(2);
            tilde_eta = para(3);
            S = para(4);
            Vdot=S'*Rpsi*(tau-M*v_r_dot-C*v_r-D*v_r+Rpsi'*Kp1*tilde_eta)-S'*Dstar*S-tilde_eta'*Kp1*Lambda*tilde_eta;
        end
        
    end
end

