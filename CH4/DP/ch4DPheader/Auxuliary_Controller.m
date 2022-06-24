classdef (Abstract) Auxuliary_Controller < Controller
    methods(Abstract)
         Vdot = lyapunov_derivative(X,U)
    end 
end