classdef nonlinearPD_controller < Auxuliary_Controller
    properties
        Kp 
        Kd
    end
    
    methods
        function obj = nonlinearPD_controller(Kp,Kd)
            obj.Kp = Kp;
            obj.Kd = Kd;
        end
        
        function U = calc_control(obj,ref,state)            
            eta_d = ref(1:3); % [xR;yR;psiR]
            eta = state(1:3); % [x;y;psi]
            vel = state(4:6); % [u;v;r]
            psi = state(3);
            R1=[cos(psi);sin(psi);0];
            R2=[-sin(psi);cos(psi);0];
            R3=[0;0;1];
            Rpsi=[R1,R2,R3];
            tilde_eta = eta - eta_d;
            eta_dot = Rpsi * vel;            
            U = -Rpsi'*(obj.Kp * tilde_eta + obj.Kd * eta_dot);
        end
                
        function get_params(obj)
            fprintf('Kp: \n')
            disp(obj.Kp)
            fprintf('Kd: \n')
            disp(obj.Kd)            
        end
        
        function Vdot = lyapunov_derivative(obj,X,tau,eta_d)
            psi = X(3);
            eta = X(1:3,1);
            vel = X(4:6,1);
            C = CV(vel);
            D = DV(vel);
            Rpsi=Rot(psi);          
            tilde_eta = eta - eta_d; 
            Vdot = vel'*(tau-C*vel-D*vel+Rpsi'*obj.Kp*tilde_eta);
        end
            
    end
end

